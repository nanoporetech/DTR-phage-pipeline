import os,sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("kmer_umap", help="File containing 2D embedding of k-mer frequency vectors (read, length, qscore, D1, D2).", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output plot file name [kmer.umap.qscore.png]", type=str, default="kmer.umap.qscore.png")
    parser.add_argument("-s", "--size", help="Size of markers [7]", type=int, default=7)
    parser.add_argument("-a", "--alpha", help="Transpancy of markers [0.2]", type=float, default=0.2)
    
    # Parse arguments
    args = parser.parse_args()

    return args

def remove_top_right_axes( ax ):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis='both', direction='in')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def scatterplot(df, title, args):
    cmap = plt.get_cmap("cool")

    fig          = plt.figure(figsize=[8,8])
 
    ax  = fig.add_axes([0.06, 0.03, 0.95, 0.8])
    res = []
    
    plot = ax.scatter(df["D1"], df["D2"], \
                      s=args.size, \
                      edgecolor="None", \
                      c=df["qscore"], \
                      cmap=cmap, \
                      alpha=args.alpha)
    
    remove_top_right_axes(ax)

    ax.set_title(title)

    fig.colorbar(plot, ax=ax, label="qscore")

    plt.xlim([df["D1"].min()-1, df["D1"].max()+1])
    plt.ylim([df["D2"].min()-1, df["D2"].max()+1])

    plt.savefig(args.output, dpi=150)

def main(args):
    df = pd.read_csv(args.kmer_umap, delimiter="\t")

    title  = "{} reads > {} bp".format(df.shape[0], min(df["length"]))

    print("Mapping {} data points...".format(df.shape[0]))
    scatterplot(df, title, args)

if __name__=="__main__":
    args = parse_args()

    main(args)