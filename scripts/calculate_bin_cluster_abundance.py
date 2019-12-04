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
    parser.add_argument('dtr_table', help='Existing statistics for each polished genome', type=str)
    parser.add_argument('bins_dir', help='Path to directory containing reads in each k-mer bin', type=str)
    parser.add_argument('clust_dir', help='Path to directory containing polished genomes for each alignment cluster', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output file containing combined CDS statistics [cds.stats.combined.tsv]', type=str, default='cds.stats.combined.tsv')
    parser.add_argument('-p', '--prefix', help=' [prefix]', type=str, default='prefix')

    # Parse arguments
    args = parser.parse_args()

    return args

def fetch_total_bin_reads(bin_id, bins_dir):
    fn      = os.path.join(bins_dir, str(bin_id), "read_list.txt")
    n_total = len(open(fn, "r").readlines())
    return n_total

def fetch_cluster_fl_reads(cluster_id, clust_dir):
    fn      = os.path.join(clust_dir, cluster_id, "{}.readinfo.csv".format(cluster_id))
    n_clust = pd.read_csv(fn).shape[0]
    return n_clust

def est_cluster_total_reads(df):
    return int(df["bin_reads"] * float(df["fl_cluster_reads"]) / df["fl_bin_reads"])

def calc_gc(cluster_id, clust_dir):
    fn           = os.path.join(clust_dir, cluster_id, "{}.pol_reads.fa".format(cluster_id))
    combined_seq = "".join([str(seq) for header,seq in SimpleFastaParser(open(fn))])
    GC           = float(combined_seq.count("G") + combined_seq.count("C")) / len(combined_seq)
    return GC

def main(args):
    dtr_df = pd.read_csv(args.dtr_table, sep="\t")
    dtr_df["sample_id"] = args.prefix
    dtr_df = dtr_df.rename(columns={"read_length":"ref_length"})

    dtr_df["bin_id"] = dtr_df.apply(lambda row: int(row["cluster_id"].split("_")[0]), axis=1)

    # Overall bin reads
    dtr_df["bin_reads"]    = dtr_df["bin_id"].apply(lambda x: fetch_total_bin_reads(x, args.bins_dir))

    # Full length bin reads
    dtr_df["fl_bin_reads"] = dtr_df["bin_id"].apply(lambda x: fetch_total_bin_reads(x, args.bins_dir))

    # Full length cluster reads
    dtr_df["fl_cluster_reads"] = dtr_df["cluster_id"].apply(lambda x: fetch_cluster_fl_reads(x, args.clust_dir))

    # Overal cluster estimated reads
    dtr_df["cluster_reads"] = dtr_df.apply(est_cluster_total_reads, axis=1)

    # Cluster reads GC content
    dtr_df["gc"] = dtr_df["cluster_id"].apply(lambda x: calc_gc(x, args.clust_dir))

    dtr_df = dtr_df.set_index("cluster_id")
    df     = dtr_df.loc[:,["bin_id","sample_id","ref_length","has_dtr","dtr_length","n_pol_reads","bin_reads","cluster_reads","gc"]]

    df.to_csv(args.output, sep="\t", index=True)

if __name__=='__main__':
    args = parse_args()

    main(args)