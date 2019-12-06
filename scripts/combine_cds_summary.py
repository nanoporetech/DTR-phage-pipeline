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
    parser.add_argument('pol_dir', help='Path to directory containing the polished genomes and annotated CDS info for each alignment cluster', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output file containing combined CDS statistics [cds.stats.combined.tsv]', type=str, default='cds.stats.combined.tsv')

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    fns = glob(os.path.join(args.pol_dir, "*", "*.cds.stats.txt"))

    df = pd.concat([pd.read_csv(fn, sep="\t").set_index("clust_id") for fn in fns], axis=0)

    df.reset_index().to_csv(args.output, sep="\t", index=False)

if __name__=='__main__':
    args = parse_args()

    main(args)