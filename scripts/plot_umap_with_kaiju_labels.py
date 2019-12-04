import os,sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from itertools import groupby
import pandas as pd
from collections import OrderedDict
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("kmer_umap", help="File containing 2D embedding of k-mer frequency vectors (read, length, qscore, D1, D2).", type=str)
    parser.add_argument("kaiju", help="Taxonomic annotation results", type=str)

    # Optional arguments
    parser.add_argument("-r", "--tax_rank", help="Taxonomic rank to plot [3]", type=int, default=3)
    parser.add_argument("-o", "--output", help="Output file [kmer.umap.kaiju.png]", type=str, default="kmer.umap.kaiju.png")
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

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

def scatterplot(df, title, rank, args):
    fig  = plt.figure(figsize=[15,8])
    ax   = fig.add_axes([0.05, 0.05, 0.5, 0.9])
    
    df   = df.sample(frac=1).reset_index(drop=True)

    # plot all the reads without labeling
    common_props = {"marker":"o", \
                    "s":args.size, \
                    "linewidth":0, \
                    "alpha":args.alpha}

    ax.scatter(df["D1"], \
               df["D2"], \
               c=df["color"], \
               label=None,
               **common_props)

    # replot one point from each lineage for the legend
    lineage_s = sorted(df["lineage"].unique())
    for i,s in enumerate(lineage_s):
        df_ = df[df["lineage"]==s].reset_index(drop=True)
        plot = ax.scatter(df_.iloc[0]["D1"], \
                          df_.iloc[0]["D2"], \
                          c=df_.iloc[0]["color"], \
                          label=s, \
                          **common_props)

    handles,labels = ax.get_legend_handles_labels()
    leg_tup = list(zip(handles, labels))
    leg_tup.sort(key=lambda x: x[1])

    # Put "unclassified" at end of legend
    unclassified = [tup for tup in leg_tup if tup[1]=="unclassified"][0]
    leg_tup.remove(unclassified)
    leg_tup.append(unclassified)
    handles, labels = zip(*leg_tup)

    ncol = 1 if len(labels)<=36 else 2
    leg = ax.legend(handles, labels, \
                    loc='upper left', \
                    bbox_to_anchor=(1, 1), \
                    prop={'size':10}, \
                    frameon=False, \
                    scatterpoints=1, \
                    ncol=ncol)

    for i in range(len(labels)):
        leg.legendHandles[i]._sizes = [100]
        leg.legendHandles[i].set_alpha(1)

    ax.set_xlim([df["D1"].min()-1, df["D1"].max()+1])
    ax.set_ylim([df["D2"].min()-1, df["D2"].max()+1])
    
    remove_top_right_axes(ax)

    ax.set_title(title)

    fig.savefig(args.output, dpi=150)

def map_annotations_to_read_ids(umap_df, labels_df):
    new_df = pd.merge(umap_df, labels_df, on="read", how="inner")
    return new_df

def add_label_details(df):
    colors  = plt.get_cmap('nipy_spectral')(np.linspace(0.03, 0.95, len(df["lineage"].unique())))
    markers = ["v", "^", ">", "<", "s", "D"]
    
    df = df.sort_values("lineage")

    for i,lin in enumerate(df["lineage"].unique()):
        df.loc[df["lineage"]==lin, "color"]  = matplotlib.colors.to_hex(colors[i], keep_alpha=True)
        df.loc[df["lineage"]==lin, "marker"] = markers[i%len(markers)]

    df.loc[df["lineage"]=="unclassified", "color"] = "gray"
    df.loc[df["lineage"]=="unclassified", "marker"] = "o"

    return df

def get_desired_lineage_rank(kaiju_annot_df, rank):
    """
    0:  cellular organisms
        Viruses
    1:  Archaea
        Bacteria
        dsDNA viruses, no RNA stage
        environmental samples
        Eukaryota
        Retro-transcribing viruses
        unclassified bacterial viruses
        unclassified viruses
    2:  LOTS
    """
    kaiju_annot_df["lineage"] = kaiju_annot_df["lineage"].replace(np.nan, "unclassified; "*20)
    kaiju_annot_df["lineage"] = kaiju_annot_df["lineage"].replace("", "unclassified; "*20)

    # Make sure there are entries at all ranks, even if classification stops at root
    new_labs = []
    for i,entry in enumerate(kaiju_annot_df["lineage"]):
        entry+="unclassified; "*20
        new_labs.append(entry)

    kaiju_annot_df["lineage"] = list(map(lambda x: x.split("; ")[rank], new_labs))
    return kaiju_annot_df

def main(args):
    umap_df        = pd.read_csv(args.kmer_umap, delimiter="\t")
    kaiju_annot_df = pd.read_csv(args.kaiju, delimiter="\t", names=["classified","read","taxid","lineage"])
    min_lab_n      = 300

    kaiju_annot_df = get_desired_lineage_rank(kaiju_annot_df, args.tax_rank)

    labels_df = kaiju_annot_df[["read","lineage"]]

    df = map_annotations_to_read_ids(umap_df, labels_df)

    # Replace any labels with insufficent counts
    df.loc[df.groupby("lineage")["lineage"].transform("count").lt(min_lab_n), "lineage"] = "unclassified"   

    # print("Total counts (N>{}):".format(min_lab_n))
    # print(df.groupby("lineage")["lineage"].agg("count"))

    df = add_label_details(df)

    title  = "{nreads} reads > {minlen} bp".format(nreads=df.shape[0], minlen=min(df["length"]))

    print("Mapping {} data points...".format(df.shape[0]))
    scatterplot(df, title, args.tax_rank, args)

if __name__=="__main__":
    args = parse_args()

    main(args)