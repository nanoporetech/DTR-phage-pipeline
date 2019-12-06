import os,sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('cds', help='Fasta file containing CDS annotated by Prodigal', type=str)
    parser.add_argument('genome', help='Fasta file containing corresponding genome sequence', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output file containing CDS stats [cds.stats.tsv]', type=str, default='cds.stats.tsv')

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    cluster_id = os.path.basename(args.cds).split('.')[0]

    cds_bases = 0
    cds_count = 0

    for header,seq in SimpleFastaParser(open(args.genome)):
        ref_len = len(seq)

    with open(args.output, 'w') as f:
        f.write('clust_id\tread\tref_len\tcds_n\tcds_sum\tcds_mean\tcds_frac\n')
        
        for header,seq in SimpleFastaParser(open(args.cds)):
            ref_read   = header.split(' ')[0].split('_')[0]
            cds_bases += len(seq)
            cds_count += 1

        cds_mean = float(cds_bases) / cds_count
        cds_frac = float(cds_bases) / ref_len
        f.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\n'.format(cluster_id,ref_read,ref_len,cds_count,cds_bases,cds_mean,cds_frac))

if __name__=='__main__':
    args = parse_args()

    main(args)