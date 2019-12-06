import os,sys
from glob import glob
import re
import numpy as np
from Bio import SeqIO
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('dtr_table', help='DTR statistics for each polished genome', type=str)
    parser.add_argument('cds_table', help='CDS statistics for each polished genome', type=str)
    parser.add_argument('bins_dir', help='Path to directory containing reads in each k-mer bin', type=str)
    parser.add_argument('clust_dir', help='Path to directory containing polished genomes for each alignment cluster', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output file containing combined CDS statistics [cds.stats.combined.tsv]', type=str, default='cds.stats.combined.tsv')
    parser.add_argument('-p', '--prefix', help=' [prefix]', type=str, default='prefix')

    # Parse arguments
    args = parser.parse_args()

    return args

def fetch_total_bin_reads(bin_id, bins_dir):
    fn      = os.path.join(bins_dir, str(bin_id), 'read_list.txt')
    n_total = len(open(fn, 'r').readlines())
    return n_total

def fetch_cluster_fl_reads(clust_id, clust_dir):
    fn      = os.path.join(clust_dir, clust_id, '{}.readinfo.csv'.format(clust_id))
    n_clust = pd.read_csv(fn).shape[0]
    return n_clust

def est_cluster_total_reads(df):
    return int(df['bin_reads'] * float(df['fl_cluster_reads']) / df['fl_bin_reads'])

def count_pol_reads(clust_id, clust_dir):
    fn = os.path.join(clust_dir, clust_id, '{}.pol_reads.fa'.format(clust_id))
    return len([Seq for Seq in SeqIO.parse(fn, 'fasta')])

def calc_gc(clust_id, clust_dir):
    fn           = os.path.join(clust_dir, clust_id, '{}.pol_reads.fa'.format(clust_id))
    combined_seq = ''.join([str(seq) for header,seq in SimpleFastaParser(open(fn))])
    GC           = float(combined_seq.count('G') + combined_seq.count('C')) / len(combined_seq)
    return GC

def main(args):
    cds_df = pd.read_csv(args.cds_table, sep='\t').drop('ref_len', axis=1).set_index('read')

    dtr_df = pd.read_csv(args.dtr_table, sep='\t').set_index('read')
    dtr_df['sample_id'] = args.prefix
    dtr_df = dtr_df.rename(columns={'seq_len':'ref_length'})

    dtr_df = dtr_df.join(cds_df, how='inner')

    dtr_df['bin_id'] = dtr_df.apply(lambda row: int(row['clust_id'].split('_')[0]), axis=1)

    # Overall bin reads
    dtr_df['bin_reads'] = dtr_df['bin_id'].apply(lambda x: fetch_total_bin_reads(x, args.bins_dir))

    # Full length cluster reads
    dtr_df['cluster_reads'] = dtr_df['clust_id'].apply(lambda x: fetch_cluster_fl_reads(x, args.clust_dir))

    # Cluster reads GC content
    dtr_df['gc'] = dtr_df['clust_id'].apply(lambda x: calc_gc(x, args.clust_dir))

    # Cluster polishing reads count
    dtr_df['n_pol_reads'] = dtr_df['clust_id'].apply(lambda x: count_pol_reads(x, args.clust_dir))

    df     = dtr_df.loc[:,['bin_id','clust_id','sample_id','ref_length','has_dtr','dtr_length','bin_reads','cluster_reads','n_pol_reads','gc']]

    df.to_csv(args.output, sep='\t', index=True)

if __name__=='__main__':
    args = parse_args()

    main(args)