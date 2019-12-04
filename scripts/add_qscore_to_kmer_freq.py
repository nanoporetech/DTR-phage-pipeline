import sys
import pandas as pd
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("kmer_freq", help="File containing k-mer frequency vectors for reads (read, length, kmer1, kmer2, etc).", type=str)
    parser.add_argument("summary", help="Sequencing summary file corresponding to reads", type=str)

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output file name [kmer_freq.qs.tsv]", type=str, default="kmer_freq.qs.tsv")

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    freq = pd.read_csv(args.kmer_freq, sep="\t")
    summ = pd.read_csv(args.summary, sep="\t", quoting=3)

    # Sometimes these summary files have quotation marks around each row entry
    if type(summ.iloc[0,0])==str:
        if summ.iloc[0,0].find("\"")>-1:
            summ.columns    = summ.columns.map(lambda x: x.strip("\""))
            summ.iloc[:,0]  = summ.iloc[:,0].str.strip("\"")
            summ.iloc[:,-1] = summ.iloc[:,-1].str.strip("\"").astype("float")

    # Adjust the column names as necessary for 1Dsq reads
    qscore = 'mean_qscore_template'
    if "sequence_length_2d" in summ.columns:
        qscore = qscore.replace('template', '2d')

    summ   = summ.rename(index=str, columns={qscore: "qscore", "read_id": "read"})

    summ = summ.set_index("read")
    freq = freq.set_index("read")

    motifs = list(freq.columns.values[2:])

    freq = freq.join(summ["qscore"], how="left").reset_index()

    all_cols = ["read", "qscore", "length"] + motifs
    freq.loc[:,all_cols].to_csv(args.output, sep="\t", index=False)

if __name__=="__main__":
    args = parse_args()

    main(args)