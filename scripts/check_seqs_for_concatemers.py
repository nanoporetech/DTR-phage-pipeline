import os,sys
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import glob
import pandas as pd
import multiprocessing
import re
import numpy as np
from tqdm import tqdm
import math
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('fasta', help='File containing the sequences to check for concatemers', type=str)

    # Optional arguments
    parser.add_argument('-o', '--output', help='TSV output file name with concatemer results [output.concats.tsv]', type=str, default='output.concats.tsv')
    parser.add_argument('-l', '--ovlp_len', help='Subsequence length to align from each end of sequence [3000]', type=int, default=3000)
    parser.add_argument('-d', '--tmp_dir', help='Temporary directory to use for writing alignments [tmp_aligns]', type=str, default='tmp_aligns')
    parser.add_argument('-t', '--threads', help='Number of threads to use [16]', type=int, default=16)
    parser.add_argument('-m', '--min_len', help='Minimum input sequence length to query for concatemers [15000]', type=int, default=15000)

    # Parse arguments
    args = parser.parse_args()

    return args


def chunks( l, n ):
    """
    Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]

def launch_pool( procs, funct, args ):
    p    = multiprocessing.Pool(processes=procs)
    try:
        results = p.map(funct, args)
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results

def run_aligner( tup ):
    orig_Seq = tup[0]
    tmp_dir  = tup[1]
    ovlp_len = tup[2]
    n        = tup[3]

    read_id   = orig_Seq.name
    # m         = re.search(r'cluster=(\d+_\d+)', orig_Seq.description)
    # if m:
    #     cluster = m.group(1)
    cluster   = '_'.join(read_id.split('_')[-2:])
    prefix    = os.path.join(tmp_dir,read_id)

    seq_start = orig_Seq.seq[:ovlp_len]
    record    = SeqRecord(seq=seq_start, id=orig_Seq.id+'.first{}'.format(ovlp_len), description=orig_Seq.id+'.first{}'.format(ovlp_len))
    start_fn  = prefix+'.first{}.fa'.format(ovlp_len)
    SeqIO.write(record, start_fn, 'fasta')

    seq_whole = orig_Seq.seq
    record    = SeqRecord(seq=seq_whole, id=orig_Seq.id, description=orig_Seq.id)
    whole_fn  = prefix+'.whole.fa'
    SeqIO.write(record, whole_fn, 'fasta')

    minimap_out = prefix+'.paf'
    minimap_CMD = 'minimap2 -x ava-ont %s %s > %s 2> /dev/null' % (start_fn, whole_fn, minimap_out)
    os.system(minimap_CMD)

    df = pd.read_csv(minimap_out, sep='\t', names=['qname','qlen','S2','E2','strand','rname','rlen','S1','E1','matches','alen1','mapqv', 'v', 'w', 'x', 'y', 'z'])
    os.remove(start_fn)
    os.remove(minimap_out)

    df = df[df['alen1']>=0.9*ovlp_len]
    df = df.sort_values('S2')

    if df.shape[0]>2:
        diffs = []
        for i,(idx,row) in enumerate(df.iterrows()):
            if i==0:
                prev_S2 = df.iloc[i,2]
                continue
            S2      = df.iloc[i,2]
            diff    = S2 - prev_S2
            prev_S2 = df.iloc[i,2]

            diffs.append(diff)
        repeat_size = np.median(diffs)
        readlen = df.iloc[0,1]
        count   = readlen/repeat_size
        return read_id,readlen,repeat_size,count

def main(args):
    func_args = []
    for i,orig_Seq in enumerate(tqdm(SeqIO.parse(args.fasta, 'fasta'), desc='Filtering reads')):
        if len(orig_Seq.seq)>=args.min_len:
            func_args.append( (orig_Seq, args.tmp_dir, args.ovlp_len, i) )

    chunk_size = int(math.ceil(len(func_args)/100))
    arg_chunks = list(chunks(func_args, chunk_size))
    results = []
    for arg_chunk in tqdm(arg_chunks, desc='Running Minimap2'):
        chunk_results = launch_pool(args.threads, run_aligner, arg_chunk)
        chunk_results = [(result[0],result[1],result[2],result[3]) for result in chunk_results if result!=None]

        results.append(chunk_results)

    with open(args.output, 'w') as f:
        f.write('read\tlength\trepeat_size\tcopies\n')
        for chunk_results in results:
            for (read_id,readlen,repeat_size,count) in chunk_results:
                f.write('{}\t{}\t{}\t{}\n'.format(read_id,readlen,repeat_size,count))

if __name__=='__main__':
    args = parse_args()

    main(args)