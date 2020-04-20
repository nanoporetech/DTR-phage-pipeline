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
    parser.add_argument('ref_read_list', help='File with list of reference reads for the  cluster', type=str)
    parser.add_argument('pol_reads_list', help='List of reads in  cluster', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output FASTA file [polished.renamed.fasta]', type=str, default='polished.renamed.fasta')

    # Parse arguments
    args = parser.parse_args()

    return args

def get_clust_info(args, clust_id):
    """ parse the read list csv file for cluster info """
    # number of polishing reads
    n_pols = len(open(args.pol_reads_list).readlines())
    # reference read
    ref_read = open(args.ref_read_list).readlines()[-1].strip()
    return ref_read, n_pols

def main(args):
    # get the cluster ID from the path name
    clust_id = os.path.basename(args.fasta).split('.')[0]

    ref_read, n_pols = get_clust_info(args, clust_id)

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
