from pcloud_clustering import pcloud_clustering 
from counting_oligos import count_oligo
from pathlib import Path
import argparse
def main():
    parser = argparse.ArgumentParser(
        prog='P-cloud clustering',
        description='Given a dictionary of oligos and their occurences, perform P-cloud clustering.'
    )

    parser.add_argument('--seq_id', type=str, required=True, help='Sequence id')
    parser.add_argument('--output', type=str, required=True, help='Path to save the pcloud results pickle file')
    parser.add_argument('--lower', type=int, default=2, help='Minimum copies for an oligo to be included in the outer layer')
    parser.add_argument('--core', type=int, default=5, help='Minimum copies for an oligo to be included in the core set')
    parser.add_argument('--secondary', type=int, default=25, help='If core oligo in P-cloud > secondary copies 2nt sufficient for inclusion in outer layer')
    parser.add_argument('--tertiary', type=int, default=50, help='If core oligo in P-cloud > tertiary copies; 3nt differences to be included in outer layer')

    args = parser.parse_args()

    cutoffs_args = {
        "lower": args.lower,
        "core": args.core,
        "secondary": args.secondary,
        "tertiary": args.tertiary
    }
    # generate the json file of k-mer counts
    dict_json_file = f'../../data/{args.seq_id}_oligo_occurences.json' 
    print(dict_json_file)
    if not Path(dict_json_file).exists():
        count_oligo(args.seq_id)

    pcloud_list = pcloud_clustering(dict_json_file, args.output, **cutoffs_args) 
    
if __name__ == "__main__":
    main()
