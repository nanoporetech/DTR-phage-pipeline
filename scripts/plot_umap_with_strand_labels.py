import os,sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.font_manager as font_manager
from matplotlib import rc,rcParams
import numpy as np
import random
from itertools import groupby
from collections import Counter
import pandas as pd
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("umap", help="UMAP embedding of k-mer counts vector for each read (read, length, qscore, D1, D2).", type=str)
    parser.add_argument("strands", help="Strand annotations gathered from each read used to polish a representative genome from an alignment cluster", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output plot file name [kmer.umap.rl.png]", type=str, default="kmer.umap.rl.png")
    
    # Parse arguments
    args = parser.parse_args()

    return args

def remove_top_right_axes( ax ):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis='both', direction='in')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def scatterplot(df, args, title):
    fig   = plt.figure(figsize=[8,8])
 
    ax  = fig.add_axes([0.06, 0.05, 0.85, 0.85])
    
    def map_colors(lab):
        if lab=="Not included":
            c = "grey"
        elif lab=="+":
            c = "r"
        elif lab=="-":
            c = "b"
        elif lab=="ref":
            c = "k"
        return c

    df = df.fillna("Not included")

    df["colors"] = df["label"].map(map_colors)

    mk = df["label"]=="Not included"
    ax.scatter(df[mk]["D1"], \
               df[mk]["D2"], \
                s=7, \
                edgecolor="None", \
                c=df[mk]["colors"], \
                alpha=0.2)

    mk = df["label"]!="Not included"
    ax.scatter(df[mk]["D1"], \
               df[mk]["D2"], \
                s=7, \
                edgecolor="None", \
                c=df[mk]["colors"], \
                alpha=0.2)

    mk = df["label"]=="ref"
    ax.scatter(df[mk]["D1"], \
               df[mk]["D2"], \
                s=3, \
                edgecolor="None", \
                c=df[mk]["colors"], \
                alpha=1.0) 

    biggest = max([max(abs(df["D1"])), max(abs(df["D2"]))])+5
    plt.xlim([-biggest,biggest])
    plt.ylim([-biggest,biggest])

    remove_top_right_axes(ax)
    ax.set_title(title)

    plt.savefig(args.output, dpi=150)

def main(args):
    umap    = pd.read_csv(args.umap, sep="\t").set_index("read")
    strands = pd.read_csv(args.strands, sep="\t").rename(columns={"read_id":"read"}).set_index("read")

    df = umap.join(strands, how="left")

    title  = "{} reads > {} bp".format(df.shape[0], min(df["length"]))

    print("Mapping {} data points...".format(df.shape[0]))
    scatterplot(df, args, title)

if __name__=="__main__":
    args = parse_args()

    main(args)