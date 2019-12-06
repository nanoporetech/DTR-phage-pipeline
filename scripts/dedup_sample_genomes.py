import os,sys
import pandas as pd
from Bio import SeqIO
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('coords', help='Input *.coords file produced by the nucmer command show-coords -c -l out.delta', type=str)
    parser.add_argument('in_fasta', help='Input fasta file containing polished genomes to deduplicate', type=str)
    parser.add_argument('in_stats', help='Input tsv file containing statistics for all polished genomes', type=str)

    # Optional arguments
    parser.add_argument('-f', '--out_fasta', help='Output fasta file name [dedup.seqs.fasta]', type=str, default='dedup.seqs.fasta')
    parser.add_argument('-s', '--out_stats', help='Output statistics file name [dedup.stats.tsv]', type=str, default='dedup.stats.tsv')
    parser.add_argument('-p', '--min_pct', help='Minimum allowable percent identity [0.98]', type=float, default=0.98)
    parser.add_argument('-c', '--min_cov', help='Minimum allowable percent genome overlap [0.98]', type=float, default=0.98)

    # Parse arguments
    args = parser.parse_args()

    return args

def remove_self_self(df, args):
    df = df.query('readr!=readq')
    df = df[(df['idy']>=args.min_pct) & \
            (df['covq']>=args.min_cov) & \
            (df['covr']>=args.min_cov)]
    return df

# def xxx():


def main(args):
    nuc_df = pd.read_csv(args.coords, skiprows=4, header=None, sep='\t')
    nuc_df = nuc_df.rename(columns={0:'S1', 1:'E1', 2:'S2', 3:'E2', 4:'len1', 5:'len2', 6:'idy', 
                                    7:'lenR', 8:'lenQ', 9:'covr', 10:'covq', 11:'readr', 12:'readq'})

    df = pd.read_csv(args.in_stats, sep='\t').set_index('read')

    # remove self-self hits and apply alignment thresholds
    nuc_df = remove_self_self(nuc_df, args)

    # Ingest all SeqIO objects
    seqs = {}
    for Seq in SeqIO.parse(args.in_fasta, 'fasta'):
        seqs[Seq.id] = Seq

    drop_read    = []
    # Figure out which ones to keep based on all vs. all aligns and Prodigal stats
    for idx,row in nuc_df.iterrows():
        # Each row represents a duplicate, pick version with most n_pol_reads
        n_pol_r = df.loc[row['readr'],'n_pol_reads']
        n_pol_q = df.loc[row['readq'],'n_pol_reads']
        
        if n_pol_q > n_pol_r:
            # Pick the query version instead of query
            drop_read.append(row['readr'])
        elif n_pol_q <= n_pol_r:
            # Pick the reference version instead of query
            drop_read.append(row['readq'])

    # Drop appropriate reads from the stats dataframe
    df = df[(~df.index.isin(drop_read)) & \
            (df.index.isin(seqs.keys()))]

    # Output the filtered set of draft genomes
    keep_fasta = []
    for seq_id,Seq in seqs.items():
        if seq_id in df.index.values:
            keep_fasta.append(Seq)

    with open(args.out_fasta, 'w') as f:
        SeqIO.write(keep_fasta, f, 'fasta')

    df.to_csv(args.out_stats, sep='\t')


if __name__=='__main__':
    args = parse_args()

    main(args)