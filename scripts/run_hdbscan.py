import os,sys
import numpy as np
import hdbscan
from sklearn import decomposition
import random
import pandas as pd
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("map_coords", help="File containing data coordinates to be binned (e.g. output from UMAP).", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output file [hdbscan.bins.tsv]", type=str, default="hdbscan.bins.tsv")
    parser.add_argument("-c", "--min_cluster", help="Minimum number of reads to call a bin [30]", type=int, default=30)

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    df = pd.read_csv(args.map_coords, delimiter="\t")
    
    X = df.loc[:,["D1", "D2"]]

    print("Running HDBSCAN to label {} UMAP points".format(df.shape[0]))
    df["bin_id"] = hdbscan.HDBSCAN(min_cluster_size=args.min_cluster).fit_predict(X)
    
    print("Created {} bins with > {} reads".format(len(df["bin_id"].unique()), args.min_cluster))
    df.to_csv(args.output, sep="\t", index=False)

if __name__=="__main__":
    args = parse_args()

    main(args)