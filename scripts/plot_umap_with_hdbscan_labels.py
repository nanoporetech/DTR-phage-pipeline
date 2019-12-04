import os,sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("hdbscan_bins", help="File containing binned data coordinates (e.g. output from HDBSCAN)", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output file [kmer.umap.hdbscan.png]", type=str, default="kmer.umap.hdbscan.png")
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

def combine_color_vals(row):
    if row["bin_id"]==-1:
        c = "grey"
    else:
        c = [row["color1"], row["color2"], row["color3"], row["color4"]]
    return c

def scatterplot(df, title, args):
    fig = plt.figure(figsize=[16,16])
    ax  = fig.add_axes([0.06, 0.03, 0.9, 0.9])
    
    bin_colors = 20
    colors = plt.get_cmap("gist_rainbow")(np.linspace(0.04, 0.96, bin_colors))
    
    for bin_id in range(-1,max(df["bin_id"])+1):
        idx = df[df["bin_id"]==bin_id].index
        df.loc[idx,"color1"] = colors[bin_id%bin_colors][0]
        df.loc[idx,"color2"] = colors[bin_id%bin_colors][1]
        df.loc[idx,"color3"] = colors[bin_id%bin_colors][2]
        df.loc[idx,"color4"] = colors[bin_id%bin_colors][3]

    df["color"] = df.apply(combine_color_vals, axis=1)
    df          = df.drop(["color1", "color2", "color3", "color4"], axis=1)

    plot = ax.scatter(df["D1"], df["D2"], \
                      s=args.size, \
                      edgecolor="None", \
                      c=df["color"], \
                      alpha=args.alpha)

    remove_top_right_axes(ax)

    ax.set_title(title)

    plt.xlim([df["D1"].min()-1, df["D1"].max()+1])
    plt.ylim([df["D2"].min()-1, df["D2"].max()+1])

    plot_no_text_fn = "{}".format('.'.join(args.output.split('.')[:-1])+'.notext.'+args.output.split('.')[-1])
    print("Saving {}".format(plot_no_text_fn))
    plt.savefig(plot_no_text_fn, dpi=300, transparent=True)
    
    with open(args.output.replace('.png','.tsv'), "w") as f:
        f.write("{}\t{}\t{}\n".format("bin", "x", "y"))
        for bin_id in range(-1,max(df["bin_id"])+1):
            sys.stderr.write("bin {} ({} reads)\n".format(bin_id,df[df["bin_id"]==bin_id].shape[0]))
            x = np.mean(df[df["bin_id"]==bin_id]["D1"])
            y = np.mean(df[df["bin_id"]==bin_id]["D2"])
            plt.text(x,y, str(bin_id), fontsize=8, horizontalalignment="center", verticalalignment="center")
            f.write("{}\t{}\t{}\n".format(bin_id, round(x,1), round(y,1)))

    remove_top_right_axes(ax)

    ax.set_title(title)

    plt.xlim([df["D1"].min()-1, df["D1"].max()+1])
    plt.ylim([df["D2"].min()-1, df["D2"].max()+1])

    print("Saving {}".format(args.output))
    plt.savefig(args.output, dpi=300)

    # Also make version without bin colors
    fig  = plt.figure(figsize=[16,16])
    ax   = fig.add_axes([0.06, 0.03, 0.9, 0.9])
    plot = ax.scatter(df["D1"], df["D2"], \
                      s=args.size, \
                      edgecolor="None", \
                      c="k", \
                      alpha=args.alpha)
    remove_top_right_axes(ax)

    ax.set_title(title)

    plt.xlim([df["D1"].min()-1, df["D1"].max()+1])
    plt.ylim([df["D2"].min()-1, df["D2"].max()+1])

    plot_fn = "{}".format('.'.join(args.output.split('.')[:-1])+'.nolabels.'+args.output.split('.')[-1])
    print("Saving {}".format(plot_fn))
    plt.savefig(plot_fn, dpi=300)

def main(args):
    df        = pd.read_csv(args.hdbscan_bins, delimiter="\t")

    title  = "{} reads and {} HDBSCAN bins".format(df.shape[0], len(df["bin_id"].unique()))

    print("Mapping {} data points...".format(df.shape[0]))
    scatterplot(df, title, args)

if __name__ == "__main__":
    args = parse_args()

    main(args)