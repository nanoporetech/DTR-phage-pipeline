import os,sys
from Bio import SeqIO
import pandas as pd
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('stats', help='Statistics for the polished genomes', type=str)
    parser.add_argument('genomes', help='FASTA file containing the polished genomes', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output file containing only the polished genomes that contain DTRs [polished.genomes.dtr.fasta]', type=str, default='polished.genomes.dtr.fasta')

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    df = pd.read_csv(args.stats, sep='\t')

    df = df.query('has_dtr==True')

    # Remove redundant ref reads (same read_id in multiple bins)
    df = df.sort_values(['read_id', 'n_pol_reads'], ascending=[True, False])
    df = df.drop_duplicates(subset='read_id', keep='first')
    df = df.set_index('cluster_id')

    dtr_seqs = []
    for Seq in SeqIO.parse(args.genomes, 'fasta'):
        trimmed_id = Seq.description.split(' ')[1]
        read_id    = trimmed_id.split('_')[0]
        clust_id   = "_".join(trimmed_id.split("_")[-2:])

        if clust_id in df.index.values:
            n_pol_reads = df.loc[clust_id, 'n_pol_reads']
            rlen        = len(Seq.seq)

            # rename sequence with additional info
            string          = '{r} cluster={c} length={l} polish_reads={n}'.format(r=read_id, \
                                                                                   c=clust_id, \
                                                                                   l=rlen, \
                                                                                   n=n_pol_reads)
            Seq.id          = string
            Seq.name        = ''
            Seq.description = ''
            dtr_seqs.append(Seq)

    SeqIO.write(dtr_seqs, args.output, 'fasta')

if __name__=='__main__':
    args = parse_args()

    main(args)