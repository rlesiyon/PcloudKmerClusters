import os
from pathlib import Path
import argparse

from pcloud_clustering import pcloud_clustering 
from counting_oligos import count_oligo
from repeat_region_annots import mapped_pcloud_to_genome
def main():
    parser = argparse.ArgumentParser(
        prog='P-cloud clustering',
        description='Given a dictionary of oligos and their occurences, perform P-cloud clustering.'
    )

    parser.add_argument('--seq_id', type=str, required=True, help='Sequence id')
    parser.add_argument('--save_dir', type=str, required=True, help='Directory to save the results')
    parser.add_argument('--lower', type=int, default=2, help='Minimum copies for an oligo to be included in the outer layer')
    parser.add_argument('--core', type=int, default=5, help='Minimum copies for an oligo to be included in the core set')
    parser.add_argument('--secondary', type=int, default=25, help='If core oligo in P-cloud > secondary copies 2nt sufficient for inclusion in outer layer')
    parser.add_argument('--tertiary', type=int, default=50, help='If core oligo in P-cloud > tertiary copies; 3nt differences to be included in outer layer')

    args = parser.parse_args()

    # Make the directory
    os.makedirs(args.save_dir, mode=0o777, exist_ok=True)

    cutoffs_args = {
        "lower": args.lower,
        "core": args.core,
        "secondary": args.secondary,
        "tertiary": args.tertiary
    }

    # count oligos
    print(f"Oligos occurences counting")
    oligos_dic = count_oligo(args.save_dir, args.seq_id)

    # construct p-cloud
    print(f"Pcloud clustering")
    pclouds_list = pcloud_clustering(args.seq_id, args.save_dir, **cutoffs_args) 

    # map p-clouds to genome coordinates
    print("Mapping p-clouds to Genome coordinates")
    mapped_pcloud_to_genome(pclouds_list, oligos_dic, args.seq_id, args.save_dir, args.core)
    
if __name__ == "__main__":
    main()
