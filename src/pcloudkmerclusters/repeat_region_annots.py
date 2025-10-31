import pickle
from pcloud_clustering import Pcloud
import json
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

def mapped_pcloud_to_genome(pcloud_list, oligos_dic, seq_id, save_dir, core):
    pcloud_mappings_stats = []
    for pcloud in pcloud_list:
        bandwidths, observed_each_bandwith, midpoint = map_pcloud(pcloud, oligos_dic) 
        pcloud_mappings_stats.append(
            tuple((bandwidths, observed_each_bandwith, midpoint))
        )
    # plot the distrubtions; 
    plot_distribution(pcloud_mappings_stats, seq_id, save_dir, core, percentile=50) 

def map_pcloud(pcloud, oligos_dic):
    return get_pcloud_distribution(pcloud, oligos_dic)

def get_pcloud_distribution(pcloud, oligos_dic):
    return pcloud_members_within_bins(pcloud, oligos_dic)
    
def pcloud_members_within_bins(pcloud, oligos_dic):
    # get core oligo
    core_oligo = list(pcloud.core_oligo.items())
    core_kmer, _ = core_oligo[0]

    # get core_oligo positions
    core_oligo_pos = oligos_dic[core_kmer][1]
    print(core_oligo_pos)

    # get pos for all oligos in the p-cloud
    beg_pos = [ oligos_dic[kmer][1][0] for kmer in pcloud.kmer_list]
    end_pos = [ oligos_dic[kmer][1][1] for kmer in pcloud.kmer_list]

    # total oligos
    total = len(beg_pos)

    # For p-clouds with 1 skip it
    if total <= 1:
        return [], [], [] 
    
    max_pos = max(end_pos)
    min_pos = min(beg_pos)
    print(max_pos)
    print(f"Total oligos: {len(beg_pos)}")
    print(f"Min pos: {min_pos}, Max pos {max_pos}")

    # get minpoint to start getting groups of oligos
    start_mid = core_oligo_pos[0]
    end_mid   = core_oligo_pos[1] 
    midpoint = (start_mid + end_mid)/2
    print(f"Core oligo: begin {start_mid}, end {end_mid} mid {midpoint}")

    bandwidth = 15
    # max_pos
    bandwidths = []
    observed_each_bandwith = []
    while end_mid < max_pos:

        # get all the oligos with a bandwith
        i = 0
        start_mid = max(midpoint - bandwidth, 0)
        end_mid = min(midpoint + bandwidth, max_pos) 
        for oligo in pcloud.kmer_list:
            if (
                oligos_dic[oligo][1][0] >= start_mid and 
                oligos_dic[oligo][1][1] <= end_mid
            ):
                i += 1
            
        print(f"Bandwith: {bandwidth}, start-mid: {start_mid} end-mid: {end_mid} found: {i}")
        if i != 0:
            bandwidths.append(tuple((start_mid, end_mid)))
            observed_each_bandwith.append(100*i/total)
        
        # double the bandwidth
        bandwidth *= 2
    return bandwidths, observed_each_bandwith, midpoint 

def plot_distribution(pcloud_mappings_stats, seq_id, save_dir, core, percentile=15): 
    plt.figure(figsize=(10, 6))
    x_coords = []
    y_values = []
    color = []

    # Get the counts
    for bandwidths, observed_each_bandwidth, _ in pcloud_mappings_stats:
        for (start, end), v in zip(bandwidths, observed_each_bandwidth):
            if v < percentile:
                for pos in [start, end]:
                    x_coords.append(pos)
                    y_values.append(v)

    # Get the bins, and counts for each oligos midpoints
    hist, bins = np.histogram(
        a = x_coords, 
        bins = 100, 
        weights = y_values
    )
    plt.plot(bins[1:], hist)

    plt.xlabel("Genome coordinates")
    plt.ylabel("No. of oligos within a range")

    # Save the figure
    plt.savefig(f"{save_dir}/pcloud_mapping_{seq_id}_{core}.png") 

def load_pk(pk_file):
    loaded_pclouds = None
    with open(pk_file, "rb") as f:
        loaded_pclouds = pickle.load(f)
    return loaded_pclouds 

def load_json(js_file):
    oligos_dic = None
    with open(js_file, "r") as f:
        oligos_dic = json.load(f)
    return oligos_dic 

if __name__ == "__main__":
    pk_file = "../../data/pcloud_results_5.pkl"
    js_file = "../../data/NZ_AP035801.1_oligo_occurences.json"
    seq_id = "NZ_AP035801.1"
    pclouds = load_pk(pk_file) 
    oligos_dic = load_json(js_file) 
    mapped_pcloud_to_genome(pclouds, oligos_dic, seq_id)


