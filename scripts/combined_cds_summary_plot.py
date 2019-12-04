import os,sys
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('stats', help='Polished genome statistics', type=str)

    # Optional arguments
    parser.add_argument('--output1', help='Output filename for plot including all polished genomes [polished.cds.stats.all.png]', type=str, default='polished.cds.stats.all.png')
    parser.add_argument('--output2', help='Output filename for plot including polished genomes with DTRs [polished.cds.stats.dtr.png]', type=str, default='polished.cds.stats.dtr.png')
    parser.add_argument('--output3', help='Output filename for plot including polished genomes with >=10 polishing reads [polished.cds.stats.n10.png]', type=str, default='polished.cds.stats.n10.png')
    parser.add_argument('--output4', help='Output filename for plot including polished genomes with >=20 polishing reads [polished.cds.stats.n20.png]', type=str, default='polished.cds.stats.n20.png')

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    figs = [args.output1, args.output2, args.output3, args.output4]
    df = pd.read_csv(args.stats, sep="\t")
    
    for i,fig_out in enumerate(figs):

        fig = plt.figure(figsize=[6,10])

        ax1 = fig.add_subplot(4,1,1)
        ax2 = fig.add_subplot(4,1,2)
        ax3 = fig.add_subplot(4,1,3)
        ax4 = fig.add_subplot(4,1,4)

        ax1.set_xlim([0,200])
        ax2.set_xlim([0,80000])
        ax3.set_xlim([0,1000])
        ax4.set_xlim([0,1.0])

        ax1.set_xlabel("# CDS")
        ax2.set_xlabel("CDS sum (bp)")
        ax3.set_xlabel("CDS avg size (bp)")
        ax4.set_xlabel("CDS fraction")

        if i==0:
            df = df
            ax1.set_title("All Prodigal CDS")
        elif i==1:
            df = df[df["has_dtr"]==True]
            ax1.set_title("Prodigal CDS for DTR genomes")
        elif i==2:
            df = df[(df["has_dtr"]==True) & (df["n_pol_reads"]>=10)]
            ax1.set_title("Prodigal CDS for DTR genomes w/ n_pols>=10")
        elif i==3:
            df = df[(df["has_dtr"]==True) & (df["n_pol_reads"]>=20)]
            ax1.set_title("Prodigal CDS for DTR genomes w/ n_pols>=20")

        mean_n,mean_sum,mean_mean,mean_frac = df.loc[:,["cds_n", "cds_sum", "cds_mean", "cds_frac"]].mean(axis=0)

        df["cds_n"].hist(bins=50, ax=ax1, color="darkblue", linewidth=0)
        df["cds_sum"].hist(bins=50, ax=ax2, color="royalblue", linewidth=0)
        df["cds_mean"].hist(bins=50, ax=ax3, color="dodgerblue", linewidth=0)
        df["cds_frac"].hist(bins=50, ax=ax4, color="lightblue", linewidth=0)

        ax1.text(0.65, 0.75, "Mean = {0:.1f}".format(mean_n), transform=ax1.transAxes)
        ax2.text(0.65, 0.75, "Mean = {0:.1f} bp".format(mean_sum), transform=ax2.transAxes)
        ax3.text(0.65, 0.75, "Mean = {0:.1f} bp".format(mean_mean), transform=ax3.transAxes)
        ax4.text(0.25, 0.75, "Mean = {0:.3f}".format(mean_frac), transform=ax4.transAxes)

        ax1.axvline(mean_n, color='r', linestyle='--')
        ax2.axvline(mean_sum, color='r', linestyle='--')
        ax3.axvline(mean_mean, color='r', linestyle='--')
        ax4.axvline(mean_frac, color='r', linestyle='--')

        fig.tight_layout()
        fig.savefig(fig_out)

if __name__=='__main__':
    args = parse_args()

    main(args)