import sys
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("stats", help="Concatemer stats file produced by check_seqs_for_concatemers.py", type=str)

    # Optional arguments
    parser.add_argument("-c", "--out_counts", help="Output plot of concatemer repeat counts [output.counts.png]", type=str, default="output.counts.png")
    parser.add_argument("-l", "--out_lengths", help="Output plot of concatemer repeat lengths [output.lengths.png]", type=str, default="output.lengths.png")

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    df = pd.read_csv(args.stats, sep="\t")
    df = df.rename(columns={"copies":"all"})
    df = df[df["all"]<30]

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    df["all"].plot.kde(bw_method=0.03, legend=False, color="k", ax=ax1)
    margin = 0.2

    colors = ["b","r","g","m","c"]
    for i,num in enumerate([4,5,6,7,8]):
        df_ = df[(df["all"]>(num-margin)) & (df["all"]<(num+margin))]
        if df_.shape[0]>1:
            df_["repeat_size"].plot.kde(bw_method=0.03, legend=False, color=colors[i], label="{} copies".format(num), ax=ax2)

        # add color shading to each interval around each integer in the copy number pdf
        ax1.axvspan((num-margin), (num+margin), alpha=0.5, color=colors[i])

    ax1.set_xlim([2,14])
    ax1.set_xticks(range(2,15))
    ax1.set_xlabel("Repeat counts")
    ax1.grid()
    ax1.set_title("Linear concatemer reads")
    plt.tight_layout()
    fig1.savefig(args.out_counts, dpi=150)

    # ax2.set_xticks(range(2,15),rotation=45)
    ax2.set_xlim([2500,15000])
    ax2.grid()
    ax1.set_title("Linear concatemer reads")
    ax2.set_xlabel("Repeat lengths")
    ax2.legend()
    plt.tight_layout()
    fig2.savefig(args.out_lengths, dpi=150)

if __name__=='__main__':
    args = parse_args()

    main(args)