import os
import argparse
from Bio import SeqIO
import pandas as pd
from glob import glob

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('fasta', help='Fasta file containing polished genome', type=str)
    parser.add_argument('align_clust_dir', help='Directory containing subdirectories of read info for each alignment cluster', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output FASTA file [polished.renamed.fasta]', type=str, default='polished.renamed.fasta')

    # Parse arguments
    args = parser.parse_args()

    return args

def get_clust_info(args, clust_id):
    clust_dir = os.path.join(args.align_clust_dir, clust_id)
    n_pols = pd.read_csv(os.path.join(clust_dir, '{}.pol_readlist.csv'.format(clust_id)), header=None).shape[0]
    ref_read = open(os.path.join(clust_dir, '{}.ref_readlist.csv'.format(clust_id))).read().strip('\n')
    return ref_read,n_pols

def main(args):
    clust_id = os.path.basename(args.fasta).split('.')[0]

    ref_read,n_pols = get_clust_info(args, clust_id)

    new_entry = []
    for Entry in SeqIO.parse(args.fasta, 'fasta'):
        Entry.id = '{r} cluster={c} length={l} polish_reads={n}'.format(r=ref_read, \
                                                                        c=clust_id, \
                                                                        l=len(Entry.seq), \
                                                                        n=n_pols)
        Entry.name = ''
        Entry.description = ''
        new_entry.append(Entry)

    SeqIO.write(new_entry, args.output, 'fasta')

if __name__=='__main__':
    args = parse_args()

    main(args)