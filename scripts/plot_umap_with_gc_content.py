import os,sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("kmer_umap", help="File containing 2D embedding of k-mer frequency vectors (read, length, qscore, D1, D2).", type=str)
    parser.add_argument("fastx", help="Fasta/fastq file of reads.", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output plot file name [kmer.umap.gc.png]", type=str, default="kmer.umap.gc.png")
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
    cmap = plt.get_cmap("magma")

    fig          = plt.figure(figsize=[8,8])
 
    ax  = fig.add_axes([0.06, 0.03, 0.95, 0.8])
    res = []
    
    plot = ax.scatter(df["D1"], df["D2"], \
                      s=args.size, \
                      edgecolor="None", \
                      c=df["GC"], \
                      cmap=cmap, \
                      alpha=args.alpha)
    
    remove_top_right_axes(ax)

    ax.set_title(title)

    fig.colorbar(plot, ax=ax, label="GC")

    plt.xlim([df["D1"].min()-1, df["D1"].max()+1])
    plt.ylim([df["D2"].min()-1, df["D2"].max()+1])

    plt.savefig(args.output, dpi=150)

def map_gc_content_to_reads(umap_df, gc_df):
    new_df = pd.merge(umap_df, gc_df, on="read", how="inner")
    return new_df

def gc_content(fastx_fn, ftype):
    gc = [] 
    if ftype=="fastq":
        for read_id, seq, qual in FastqGeneralIterator(open(fastx_fn)):
            seq_l   = list(seq)
            read_gc = (seq_l.count("G") + seq_l.count("C")) / float(len(seq_l))
            # read_gc = 0.5
            gc.append( (read_id.split(" ")[0], read_gc) )
    elif ftype=="fasta":
        for read_id, seq in SimpleFastaParser(open(fastx_fn)):
            seq_l   = list(seq)
            read_gc = (seq_l.count("G") + seq_l.count("C")) / float(len(seq_l))
            # read_gc = 0.5
            gc.append( (read_id.split(" ")[0], read_gc) )
            
    gc_df = pd.DataFrame(gc, columns=["read","GC"])
    return gc_df

def check_input_format(fastx_fn):
    for line in open(fastx_fn).readlines():
        break

    if line[0]=="@":
        ftype = "fastq"
    elif line[0]==">":
        ftype = "fasta"
    else:
        raise("Unexpected file type! Only *.fasta, *.fa, *.fsa, *.fna, *.fastq, and *.fq recognized.")
    return ftype

def main(args):
    umap_df  = pd.read_csv(args.kmer_umap, delimiter="\t")

    ftype = check_input_format(args.fastx)

    gc_df = gc_content(args.fastx, ftype)

    df = map_gc_content_to_reads(umap_df, gc_df)

    title  = "{} reads > {} bp".format(umap_df.shape[0], min(umap_df["length"]))

    print("Mapping {} data points...".format(df.shape[0]))
    scatterplot(df, title, args)

if __name__=="__main__":
    args = parse_args()

    main(args)