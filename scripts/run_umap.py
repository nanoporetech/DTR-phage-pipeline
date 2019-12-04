import os,sys
import numpy as np
import umap
from sklearn import decomposition
import random
import pandas as pd
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("kmer_freq", help="File containing log-normalized k-mer frequency vectors for reads (read, length, kmer1, kmer2, etc).", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output file [umap.tsv]", type=str, default="umap.tsv")
    parser.add_argument("-l", "--min_length", help="Minimum read length to include in the 2D map [15000]", type=int, default=15000)
    parser.add_argument("-q", "--min_q", help="Minimum Q-score to include in the 2D map [8]", type=int, default=8)
    parser.add_argument("-d", "--min_dist", help="Minimum distance apart that points are allowed to be in the 2D map [0.1]", type=float, default=0.1)
    parser.add_argument("-n", "--n_neighbors", help="Number of neighbors to look at when learning the manifold structure [15]", type=int, default=15)

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    df = pd.read_csv(args.kmer_freq, delimiter="\t")
    
    # Filter out short and low qscore reads
    df = df[(df["length"]>=args.min_length) & (df["qscore"]>=args.min_q)].reset_index(drop=True)

    motifs = [x for x in df.columns.values if x not in ["read", "length"]]
    X = df.loc[:,motifs]

    print("Running UMAP to embed {} features into 2 dimensions".format(X.shape[1]))
    X_embedded = umap.UMAP(n_neighbors=args.n_neighbors, min_dist=args.min_dist, verbose=2).fit_transform(X)

    df_umap = pd.DataFrame(X_embedded, columns=["D1", "D2"])

    df_to_write = pd.concat([df["read"], df["length"], df["qscore"], df_umap], axis=1)
    
    df_to_write.to_csv(args.output, sep="\t", index=False)

if __name__=="__main__":
    args = parse_args()

    main(args)