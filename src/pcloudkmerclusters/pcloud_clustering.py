'''
Assumption: 
- Equal base frequencies
- Poisson distribution

- highest frequency oligo to initiate a cloud
- Expand the cloud by adding similar high-frequency ologos to form a P-cloud `core`.
- Similar: mean differences of up to 3 nt from previosuly identified core oligo. 
    the definitions of “similar” and “high frequency” were free parameters or adjustable “cutoffs”
    - the definitions of “similar” and “high frequency” were free parameters or adjustable “cutoffs”

Pseudocode:

Initialize first P-cloud
- find the highest frequency oligo
- Add similar high-frequency oligos (copies > core cutoff) to form a P-cloud core (simillar means differences of up to 3nt from previosuly identified core oligo)

Multiple p-clouds: (core cutoff)
while (exist oligos with copies > core cutoff):
- Remove oligos belonging to identified P-cloud
- Repeat core-identification and core expansion process with the remainder (until no oligos remained with counts greater than the `core cutoff`)

Outer layer of each P-cloud
- Get all medium-copy oligos with copies == lower cutoff
- Check where it belongs as follows
- Similar varied amongst clouds: 
    - Most P-clouds a candidate oligo for the outer layer needed to have only a single difference from a core oligo. 
    - if core oligo in P-cloud > 200 copies 2nt sufficient for inclusion. (secondary cutoff)
    - if core oligo in P-cloud > 2000 copies; 3nt differences to be included. (tertiary cutoff)
    - If oligo belogs to two or more different P-clouds (it was assigned to the P-cloud with the highest frequency oligo in its core)
'''
import json
from rapidfuzz.distance import Hamming
import numpy as np
import pickle
import argparse

class Pcloud:
    '''
    P-cloud object
    core_oligo: dict of kmer and its frequency
    kmer_list: list of all kmers in the pcloud (core, core_set, outer_layer)
    core_set: dict of kmers and their frequency in the core set
    outer_layer: dict of kmers and their frequency in the outer layer
    '''
    def __init__(self, kmer, freq):
        self.core_oligo = { kmer : freq} 
        self.kmer_list = [kmer]
        self.core_set = {}
        self.outer_layer = {}
    def update_kmer_list(self):
        self.kmer_list.extend(self.core_set.keys())

def pcloud_clustering(seq_id, save_dir, lower, core, secondary, tertiary):
    '''
    Perform P-cloud clustering on the given oligos dictionary.
    - Initialize first P-cloud
    - Create multiple p-clouds until no oligos with copies > core cutoff
    - Outer layer of each P-cloud
    0. oligos_dic: dictionary of oligos and their occurences
    1. lower: minimum copies for an oligo to be included in the outer layer
    2. core: minimum copies for an oligo to be included in the core set
    3. primary: not used in the current implementation
    4. secondary: if core oligo in P-cloud > secondary copies 2nt sufficient for inclusion in outer layer
    5. tertiary: if core oligo in P-cloud > tertiary copies; 3nt differences to be included in outer layer
    return: list of Pcloud objects
    '''
    # Make the output, and dict_json_file
    dict_json_file = f"{save_dir}/{seq_id}.json" 
    output_file = f"{save_dir}/{seq_id}.pkl"

    # load oligos occurences dictionary
    with open(dict_json_file, "r") as f:
        oligos_dic = json.load(f)

    print(f'Total oligos counted: {len(oligos_dic)}')

    # Make the initial p-cloud
    pcloud = initialize_pcloud(oligos_dic, core)
    print(f'Initial P-cloud core oligo: {pcloud.core_oligo}, core set size: {len(pcloud.core_set)}')
    pcloud_list, oligos_dic = multiple_clouds(pcloud, oligos_dic, core)
    print(f'Total P-clouds created: {len(pcloud_list)}')
    pcloud_list, outer_candidates = get_outer_layer(pcloud_list, oligos_dic, lower, secondary, tertiary)
    print(f'Remaining oligo with no clusters: {len(outer_candidates)}')

    # save the pcloud results
    with open(output_file, "wb") as f:
        pickle.dump(pcloud_list, f) 
    
    return pcloud_list

def initialize_pcloud(oligos_occurences_dic, core_cutoff):
    '''
    Initialize the P-cloud by finding the highest frequency oligo amongst the oligos dictionary
    oligos_occurences_dic: dictionary of oligos and their occurences
    core_cutoff: minimum copies for an oligo to be included in the core set 
    return: Pcloud object
    '''
    # get the highest frequency oligos
    kmer = max(oligos_occurences_dic, key=lambda k: oligos_occurences_dic[k])
    freq, _ = oligos_occurences_dic[kmer]

    # constructr the Pcloud
    pcloud = Pcloud(kmer, freq)
    pcloud = expand_pcloud_core_set(pcloud, oligos_occurences_dic, core_cutoff)
    pcloud.update_kmer_list() 

    return pcloud

def expand_pcloud_core_set(pcloud, oligos_occurences_dic, core_cutoff): 
    '''
    Expand the core set of a pcloud by adding similar high-frequency ologos to form a P-cloud `core seta`.
    Similar: mean differences of up to 3 nt from previosuly identified core oligo; 
    high-frequency: copies > core cutoff
    '''
    kmers_list = [k for k, v in oligos_occurences_dic.items() if k not in pcloud.core_oligo.keys()]
    # add the core set.
    for kmer in kmers_list:
        # calculate hamming distance 
        distance = get_hamming_distance(list(pcloud.core_oligo.keys())[0], kmer)
        if distance <= 3 and oligos_occurences_dic[kmer][0] >= core_cutoff:
            pcloud.core_set[kmer] = oligos_occurences_dic[kmer]
    return pcloud

def multiple_clouds(pcloud, oligos_dic, core):
    '''
    Create multiple p-clouds until no oligos with copies > core cutoff
    '''
    # Initialize p cloud list
    pcloud_list = [pcloud]

    # Remaining oligos not in 1st pcloud
    remainder_oligos, oligos_copies, oligos_dic = update_remainder_kmer_list(pcloud, oligos_dic.keys(), oligos_dic)

    # loop through the list and make pclouds
    while (oligos_copies > core).any():
        pcloud = initialize_pcloud(oligos_dic, core)
        pcloud_list.append(pcloud)
        # update remainder_oligos, and oligos copy
        remainder_oligos, oligos_copies, oligos_dic = update_remainder_kmer_list(pcloud, remainder_oligos, oligos_dic)
        print(f'New P-cloud core oligo: {pcloud.core_oligo}, core set size: {len(pcloud.core_set)}, Remaining oligos: {len(remainder_oligos)}')

    return pcloud_list, oligos_dic

def get_outer_layer(pcloud_list, oligos_dic, lower, secondary, tertiary):
    '''
     Outer layer of each P-cloud
    - Get all medium-copy oligos with copies == lower cutoff
    - Check where it belongs as follows
    - Similar varied amongst clouds: 
        - Most P-clouds a candidate oligo for the outer layer needed to have only a single difference from a core oligo. 
        - if core oligo in P-cloud > 200 copies 2nt sufficient for inclusion. (secondary cutoff)
        - if core oligo in P-cloud > 2000 copies; 3nt differences to be included. (tertiary cutoff)
        - If oligo belogs to two or more different P-clouds (it was assigned to the P-cloud with the highest frequency oligo in its core)
    '''
    outer_candidates = [k for k, v in oligos_dic.items() if v[0] >= lower]

    # find the potential possible outer layer for the candidates.
    candidate_outer_layers = {} 
    for oligo_candidate in outer_candidates:
        for idx, pcloud in enumerate(pcloud_list):
            core_oligo = list(pcloud.core_oligo.keys())[0]
            max_freq_oligo = pcloud.core_oligo[core_oligo]


            core_oligo_distance = get_hamming_distance(oligo_candidate, core_oligo) 
            distances = [core_oligo_distance] + [get_hamming_distance(oligo_candidate, k) for k in pcloud.core_set.keys()]
            min_distance = min(distances)
            # For most P-clouds, a candidate oligo for the outer layer needed to have only a single difference from a core oligo
            if min_distance == 1: 
                if candidate_outer_layers.get(oligo_candidate):
                    candidate_outer_layers[oligo_candidate].append((pcloud, idx))
                else:
                    candidate_outer_layers[oligo_candidate] = [(pcloud, idx)]
                continue
            # If a core oligo in the P-cloud had more than 200 copies (the secondary cutoff), for example, a difference of 2 nt was sufficient for inclusion
            elif min_distance == 2 and max_freq_oligo >= secondary:
                if candidate_outer_layers.get(oligo_candidate):
                    candidate_outer_layers[oligo_candidate].append((pcloud, idx))
                else:
                    candidate_outer_layers[oligo_candidate] = [(pcloud, idx)]
            # There must be an oligo with 2000 copies in the core layer) that would allow oligos with up to 3 nt difference to be included
            elif min_distance == 3 and max_freq_oligo >= tertiary:
                if candidate_outer_layers.get(oligo_candidate):
                    candidate_outer_layers[oligo_candidate].append((pcloud, idx))
                else:
                    candidate_outer_layers[oligo_candidate] = [(pcloud, idx)]
    # Assign the candidate to their respective pclouds
    for oligo, pclouds in candidate_outer_layers.items():
        if len(pclouds) == 1:
            pcloud, idx = pclouds[0]
            pcloud.outer_layer[oligo] = oligos_dic[oligo]
        else:
            max_pcloud = max(pclouds, key=lambda x: list(x[0].core_oligo.values())[0])[0]
            max_pcloud.outer_layer[oligo] = oligos_dic[oligo]
    return pcloud_list, outer_candidates

def update_remainder_kmer_list(pcloud, remainder_list, oligos_dic):
    '''
    Update the remainder k-mer list by removing k-mers that are in the P-cloud
    '''
    remainder_oligos = [ k for k in remainder_list if k not in pcloud.kmer_list]
    oligos_copies = np.array([oligos_dic[k][0] for k in remainder_oligos])
    new_oligo_dec = { k: oligos_dic[k] for k in remainder_oligos }
    return remainder_oligos, oligos_copies, new_oligo_dec

def get_hamming_distance(kmer1, kmer2):
    '''
    Get hamming distance between two kmers
    '''
    distance = min(
        Hamming.distance(kmer1, kmer2),
        Hamming.distance(kmer1, kmer2[::-1]), 
        Hamming.distance(kmer1[::-1], kmer2)
    )
    return distance

def simulate_pclouds():
  pass


