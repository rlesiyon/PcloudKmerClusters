# Library
import math
from bitarray import bitarray
import numpy as np
from Bio import SeqIO
import os
import time
import sys
import json
from tqdm import tqdm
from pathlib import Path

# constant; 
NUCLEOTIDER_INDEX = {
    'A': 0, 
    'T': 1, 
    'G': 2,
    'C': 3
}

def count_oligo_occurences(seq_fp, method = "mixed"):
    '''
    Count the occurrences of oligonucleotides in a genomic sequence.
    Args:
        seq_fp (str): File path to the FASTA file containing the genomic sequence.
        method (str): Counting method to use. Options are "mixed_count" or "direct_count".
    Returns:
        dict or np.ndarray: A dictionary with oligonucleotide sequences as keys and their counts as values (for mixed_count),
                            or a numpy array with counts indexed by oligonucleotide indices (for direct_count).
    '''
    print('Reading genome sequence...')
    genome_sequence = read_dna_sequence(seq_fp)
    for sequence in genome_sequence: 
        n = len(sequence)
        oligo_length = get_oligo_length(n)
        if oligo_length > 16:
            print(f'Oligo length {oligo_length} is too long, setting to 16')
            oligo_length = 16
        if method == "mixed_count":
            print('Using mixed counting method...')
            return mixed_count(sequence, n, oligo_length, seq_fp)
        if method == "direct_count":
            print('Using direct counting method...')
            return direct_count(sequence, n, oligo_length, seq_fp)

def mixed_count(genome_sequence, n, oligo_length, fp): 
    '''
    Count the occurrences of oligonucleotides in a genomic sequence using a mixed counting method.
    The method useds a combination of a bit array and a hash set to efficiently count oligonucleotide occurrences.
    It saves the results to a JSON file.
    Args:
        genome_sequence (str): The genomic sequence.
        n (int): Length of the genomic sequence.
        oligo_length (int): Length of the oligonucleotides to count.
        fp (str): File path to the FASTA file (used for naming output files).
    Returns:
        dict: A dictionary with oligonucleotide sequences as keys and their counts as values.
    '''
    # initialize the counting bit array; 
    print(f'Initializing mixed counting method...')
    print(f'Oligo length: {4**oligo_length}')
    counting_bit_array = bitarray(4**oligo_length)
    print(f'Bit array size: {round(sys.getsizeof(counting_bit_array)/1024/1024/1024, 5)} Gb')
    # initialize a hash set
    oligo_occurence_hash_set = {}
    for i in tqdm(range(n - oligo_length + 1)):
        oligo_seq = genome_sequence[i:i + oligo_length]
        # keep track of the word where the index is coming from. 
        index, rep_oligo_seq = get_index_and_word(oligo_seq)
        if counting_bit_array[index] == 1: 
           if oligo_occurence_hash_set.get(str(rep_oligo_seq)):
              oligo_occurence_hash_set[str(rep_oligo_seq)][0] += 1
           else:  
              # the count is set to 2; adding the previously seen oligo, and the current one. 
              oligo_occurence_hash_set[str(rep_oligo_seq)] = [2, (i, i + oligo_length)]
        else: 
            counting_bit_array[index] = 1

    # print memory used. 
    total_size = sys.getsizeof(counting_bit_array) + sys.getsizeof(oligo_occurence_hash_set)
    print(f'Mixed method memory usage: {round(total_size/1024/1024/1024, 5)} Gb')
    
    # write to json file; 
    fp = fp.split('.fasta')[0]
    with open(f'{fp}_oligo_occurences.json', 'w') as f:
        json.dump(oligo_occurence_hash_set, f)
    return f'{fp}_oligo_occurences.json' 

def direct_count(genome_sequence, n, oligo_length, fp):
    '''
    Count the occurrences of oligonucleotides in a genomic sequence using a direct counting method.
    The method uses a numpy array to store counts of each oligonucleotide.
    Args:
        genome_sequence (str): The genomic sequence.
        n (int): Length of the genomic sequence.
        oligo_length (int): Length of the oligonucleotides to count.
        fp (str): File path to the FASTA file (not used in this method).
    Returns:
        np.ndarray: A numpy array with counts indexed by oligonucleotide indices.'''
    # initialize a counting array
    counting_array =  np.zeros(4**oligo_length)
    for i in range(n - oligo_length + 1):
        oligo_seq = genome_sequence[i:i + oligo_length]
        index, _ = get_index_and_word(oligo_seq)
        counting_array[index] += 1
    total_size = sys.getsizeof(counting_array)
    print(f'Direct method memory usage: {round(total_size/1024/1024/1024, 4)} Gb')
    return counting_array 

def read_dna_sequence(file_name):
    '''
    Read a DNA sequence from a FASTA file using SeqIO from Biopython.
    Args:
        file_name (str): Path to the FASTA file.
    '''
    with open(file_name) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield record.seq

def download_fasta_sequence(seq_id, filename):
    '''
    Download a FASTA sequence from NCBI using efetch and esearch.
    Args:
        seq_id (str): The sequence ID to download.
    ''' 
    #os.makedirs(filename, exist_ok=True)
    cmd = f'esearch -db nucleotide -query "{seq_id}" | efetch -format fasta > {filename}'
    os.system(cmd)

def get_array_index(oligo_seq): 
    '''
    Convert an oligonucleotide sequence to its corresponding index in a counting array.
    Args:
        oligo_seq (str): The oligonucleotide sequence.
    Returns:
        int: The index corresponding to the oligonucleotide sequence.
    '''
    index = 0
    for i, base in enumerate(oligo_seq):
        if base not in NUCLEOTIDER_INDEX:
            continue
        index += NUCLEOTIDER_INDEX[base]*(4**i)
    return index

def get_oligo_length(n):
    return round(math.log(n) + 1)

def get_index_and_word(oligo_seq):
    '''
    Get the index and the representative oligonucleotide sequence (smaller between the sequence and its reverse).
    Args:
        oligo_seq (str): The oligonucleotide sequence.
    Returns:
        tuple: A tuple containing the index and the representative oligonucleotide sequence.
    '''
    index1, index2 = get_array_index(oligo_seq), get_array_index(oligo_seq[::-1]) 
    if index1 < index2: 
        return index1, oligo_seq
    return index2, oligo_seq[::-1] 


def count_oligo(seq_id, method='mixed_count'):
    seq_fp = f'../../data/{seq_id}.fasta'
    if not Path(seq_fp).exists():
        download_fasta_sequence(seq_id, filename=seq_fp) 
    count_oligo_occurences(seq_fp, method = method)  

if __name__ == "__main__":
    seq_id = "AP028058.1" # Candidatus Nasuia deltocephalinicola NCIN DNA, complete genome
    #seq_id = "NZ_AP035801.1" # Escherichia coli strain E. coli JCM:5491 chromosome, complete genome
    #seq_id = "NC_000001.11" # Homo sapiens chromosome 1, complete genome
    print('Counting using mixed method...')
    start_time = time.time()
    oligo_occurence_hash_set = count_oligo(seq_id, method = "mixed_count") 
    print(f'Mixed method: {(time.time() - start_time)}')

