import os,sys
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import glob
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('fasta', help='FASTA file of polished genomes', type=str)
    parser.add_argument('pol_dir', help='Directory containing the polished genomes from each alignment cluster', type=str)

    # Optional arguments
    parser.add_argument('--overlap', help='Fraction of genome at start/end to check for DTR [0.2]', type=float, default=0.2)
    parser.add_argument('-o', '--output', help='Output file of DTR statistics [output.dtr.stats.tsv]', type=str, default='output.dtr.stats.tsv')

    # Parse arguments
    args = parser.parse_args()

    return args

def main(args):
    fasta_dir = "/".join(args.fasta.split("/")[:-1])
    pol_fas   = glob.glob(os.path.join(args.pol_dir, "*/*.pol_reads.fa"))
    pol_fa_d  = dict([ (fn.split("/")[-2], fn) for fn in pol_fas])

    # iterate over all reads in the cluster
    records = []

    for orig_Seq in SeqIO.parse(args.fasta, "fasta"):
        trimmed_id = orig_Seq.description.split(' ')[1]
        read_id   = trimmed_id.split('_')[0]
        prefix    = os.path.join(fasta_dir,read_id)
        read_len  = len(orig_Seq.seq)

        overlap_len = int(args.overlap*read_len)

        clust_id  = "_".join(trimmed_id.split("_")[-2:])
        pol_fa    = pol_fa_d[clust_id]
        n_pols    = 0
        
        for read,seq in SimpleFastaParser(open(pol_fa)):
            n_pols+=1

        seq_start = orig_Seq.seq[:overlap_len]
        record    = SeqRecord(seq=seq_start, id=orig_Seq.id+".first{}".format(overlap_len), description=orig_Seq.id+".first{}".format(overlap_len))
        start_fn  = prefix+".first{}.fa".format(overlap_len)
        SeqIO.write(record, start_fn, "fasta")

        seq_end = orig_Seq.seq[-overlap_len:]
        record  = SeqRecord(seq=seq_end, id=orig_Seq.id+".last{}".format(overlap_len), description=orig_Seq.id+".last{}".format(overlap_len))
        end_fn  = prefix+".last{}.fa".format(overlap_len)
        SeqIO.write(record, end_fn, "fasta")

        nucmer_out = "%s.out" % orig_Seq.id.split("|")[0]
        nucmer_CMD = "nucmer -p %s -c 10 %s %s > /dev/null 2>&1" % (prefix, start_fn, end_fn)
        os.system(nucmer_CMD)
        delta      = "{}.delta".format(prefix)

        coord      = "{}.coords".format(prefix)
        sc_CMD     = "show-coords -l {} > {}".format(delta, coord)
        os.system(sc_CMD)
        
        try:
            # get the sum of self-self alignments and see if those alignment
            # extend all the way to the beginning and end of the read sequence.
            # This excludes some reads that have long DTRs but also seem to be
            # concatamers (extra sequence appended at the end).
            df        = pd.read_csv(coord, skiprows=5, header=None, sep=r"\s+")
            dtr       = True
            aln_sum   = df.iloc[:,6].sum()
            aln_start = min(df.iloc[:,0])
            aln_end   = max(df.iloc[:,4])

            # transform the aln_end to reflect num. bases before end of genome/read
            aln_end   = aln_end-overlap_len

            # If the alignments are too far from the end of the read, don't call it a DTR
            if aln_start>200 or aln_end<-200:
                dtr = False

        except:
            dtr       = False
            aln_sum   = 0
            aln_start = -1
            aln_end   = 1

        os.remove(start_fn)
        os.remove(end_fn)
        os.remove(delta)
        os.remove(coord)

        records.append( (read_id, clust_id, read_len, dtr, aln_sum, aln_start, aln_end, n_pols) )

    df = pd.DataFrame.from_records(records, columns=['read_id', 'cluster_id', 'read_length', \
                                                     'has_dtr', 'dtr_length', 'dtr_start', \
                                                     'dtr_end', 'n_pol_reads'])
    df.to_csv(args.output, sep='\t', index=False)

if __name__=='__main__':
    args = parse_args()

    main(args)