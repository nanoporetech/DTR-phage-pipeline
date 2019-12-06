import sys
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("stats", help="Concatemer stats file produced by check_seqs_for_concatemers.py", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output plot of concatemer repeat lengths vs. repeat counts [output.contours.png]", type=str, default="output.contours.png")
    parser.add_argument("--min_length", help="Minimum repeat length bound [2000]", type=int, default=2000)
    parser.add_argument("--max_length", help="Maximum repeat length bound [20000]", type=int, default=20000)
    parser.add_argument("--min_copies", help="Minimum repeat count bound [3]", type=int, default=3)
    parser.add_argument("--max_copies", help="Maximum repeat count bound [8]", type=int, default=8)

    # Parse arguments
    args = parser.parse_args()

    return args

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def main(args):
    df = pd.read_csv(args.stats, sep="\t")

    sns.set(style="ticks")

    ref_sizes = [40000,67000]

    df = df[(df["copies"]>=args.min_copies) & \
            (df["copies"]<=args.max_copies) & \
            (df["repeat_size"]>=args.min_length) & \
            (df["repeat_size"]<=args.max_length)]

    fig = plt.figure(figsize=[12,12])

    X = df["repeat_size"]
    Y = df["copies"]

    orig_cmap = matplotlib.cm.Greys

    shifted_cmap = shiftedColorMap(orig_cmap, midpoint=0.30, stop=0.8, name='shifted')
    g = sns.JointGrid(X, Y)
    g = g.plot_joint(sns.kdeplot, \
                     cmap=shifted_cmap, \
                     shade=True, \
                     shade_lowest=True, \
                     n_levels=30)
    sns.kdeplot(X, color="b", shade=True, bw=0.1, vertical=False, ax=g.ax_marg_x, legend=False)
    sns.kdeplot(Y, color="b", shade=True, bw=0.1, vertical=True,  ax=g.ax_marg_y, legend=False)

    ax = g.ax_joint

    # Plot repeat size and copy number relationship for genome of <phage_size>
    x  = np.linspace(args.min_length,args.max_length,100) # 100 linearly spaced numbers
    ls = ["--", ":"]
    cs = ["k", "k"]
    for i,phage_size in enumerate(ref_sizes):
        y = phage_size/x
        ax.plot(x,y, color=cs[i], linestyle=ls[i], label="{}".format(int(phage_size/1000)))

    ax.scatter(X,Y,s=1,color="k",alpha=0.2, label='_nolegend_')

    ax.set_xlim([args.min_length,args.max_length])
    ax.set_ylim([args.min_copies,args.max_copies])
    ax.legend(loc="upper right", frameon=False, title="Sequence size (kbp)")

    ax.set_xticklabels(list(map(lambda x: int(x/1000), ax.get_xticks())))
    ax.set_xlabel("Repeat length (kbp)")
    ax.set_ylabel("Copies")

    plt.tight_layout()
    plt.savefig(args.output, dpi=200)

if __name__=='__main__':
    args = parse_args()

    main(args)