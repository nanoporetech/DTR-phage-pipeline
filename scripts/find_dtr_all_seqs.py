import os,sys
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import argparse
from tqdm import tqdm
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shutil
import multiprocessing
from collections import defaultdict

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('fastx', help='Fasta/fastq file containing read sequences', type=str)

    # Optional arguments
    parser.add_argument('-p', '--prefix', help='Output file prefix (for <prefix>.dtr.stats.tsv and <prefix>.dtr.fasta) [output]', type=str, default='output')
    parser.add_argument('-o', '--overlap', help='Check for overlaps between the first and last <overlap> percent of the sequence [20]', type=int, default=20)
    parser.add_argument('-d', '--tmpdir', help='Path to tmp directory where alignment data is written [./aln_tmp]', type=str, default='aln_tmp')
    parser.add_argument('-t', '--threads', help='Number of threads to use [4]', type=int, default=4)
    parser.add_argument('-c', '--chunk_size', help='Number of reads per chunk [5000]', type=int, default=5000)
    parser.add_argument('-m', '--hist_max', help='Maximum read length to include in histogram of DTR read lengths [80000]', type=int, default=80000)
    parser.add_argument('--no_fasta', help='Do not write FASTA of DTR sequences [False]', action='store_true')
    parser.add_argument('--no_hist', help='Do not create histogram of DTR sequence lengths [False]', action='store_true')

    # Parse arguments
    args = parser.parse_args()

    return args

def launch_pool( procs, funct, args ):
    p    = multiprocessing.Pool(processes=procs)
    try:
        results = p.map(funct, args)
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results

def chunks( l, n ):
    """
    Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]

def subseq_read_and_run_minimap2(read_id, seq, args):
    prefix    = os.path.join(args.tmpdir, read_id)
    read_len  = len(seq)

    ovlp_len = int((args.overlap/100) * read_len)

    seq_start = Seq.Seq(seq[:ovlp_len])
    record    = SeqRecord(seq=seq_start, id=read_id+'.first{}'.format(ovlp_len), description=read_id+'.first{}'.format(ovlp_len))
    start_fn  = prefix+'.first{}.fa'.format(ovlp_len)
    SeqIO.write(record, start_fn, 'fasta')

    seq_end   = Seq.Seq(seq[-ovlp_len:])
    record    = SeqRecord(seq=seq_end, id=read_id+'.last{}'.format(ovlp_len), description=read_id+'.last{}'.format(ovlp_len))
    end_fn    = prefix+'.last{}.fa'.format(ovlp_len)
    SeqIO.write(record, end_fn, 'fasta')

    minimap2_paf = '{}/{}.paf'.format(args.tmpdir, read_id.split('|')[0])
    minimap2_CMD = 'minimap2 -t 1 -x map-ont {} {} > {} 2> /dev/null'.format(end_fn, start_fn, minimap2_paf)
    os.system(minimap2_CMD)

    os.remove(start_fn)
    os.remove(end_fn)

    return minimap2_paf,ovlp_len

def write_dtr_stats(df, args):
    fname = '{}.dtr.stats.tsv'.format(args.prefix)
    print('Writing {}'.format(fname))
    df.to_csv(fname, sep='\t', index=False)

def write_dtr_fasta(df, args, ftype):
    fname = '{}.dtr.fasta'.format(args.prefix)
    print('Writing {}'.format(fname))
    dtr_ids = df[df['has_dtr']==True]['read']
    dtr_seqs = [s for s in SeqIO.parse(args.fastx, ftype) if s.id in dtr_ids.values]
    SeqIO.write(dtr_seqs, fname, 'fasta')

def check_input_format(fastx):
    for line in open(fastx).readlines():
        break

    if line[0]=="@":
        ftype = "fastq"
    elif line[0]==">":
        ftype = "fasta"
    else:
        raise("Unexpected file type! Only *.fasta, *.fa, *.fsa, *.fna, *.fastq, and *.fq recognized.")
    return ftype

def get_n_reads(fastx, ftype):
    n_lines = 0
    with open(fastx) as f:
        for i, l in enumerate(f):
            n_lines += 1
    
    if ftype=="fastq":
        n_reads = len([read_tup for read_tup in FastqGeneralIterator(open(fastx))])
    elif ftype=="fasta":
        n_reads = len([read_tup for read_tup in SimpleFastaParser(open(fastx))])
    return n_reads

def search_for_dtr( tup ):
    read_id        = tup[0]
    seq            = tup[1]
    
    paf,ovlp_len = subseq_read_and_run_minimap2(read_id, seq, args)
    df = pd.read_csv(paf, sep='\t', usecols=range(12), \
                     names = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', \
                     'tlen', 'tstart', 'tend', 'matches', 'alnlen', 'mapqv'])
    
    aln_len = 'alnlen'
    q_start = 'qstart'
    t_end   = 'tend'

    os.remove(paf)

    if df.shape[0]>0: # an alignment exists
        # get the sum of self-self alignments and see if those alignment
        # extend all the way to the beginning and end of the sequence.
        dtr       = True
        aln_sum   = df.loc[:,aln_len].sum()
        aln_start = min(df.loc[:,q_start])
        aln_end   = max(df.loc[:,t_end])

        # transform the aln_end value to reflect num. bases before end of genome/read
        aln_end   = aln_end - ovlp_len

        # If the alignments are too far from the end of the read, don't call it a DTR
        if aln_start>200 or aln_end<-200:
            dtr = False
    else: # no alignment exists
        dtr       = False
        aln_sum   = 0
        aln_start = -1
        aln_end   = 1

    read_dtr_stats = {}

    read_dtr_stats['read'] = read_id
    read_dtr_stats['seq_len'] = len(seq)
    read_dtr_stats['has_dtr'] = dtr
    read_dtr_stats['dtr_length'] = aln_sum
    read_dtr_stats['dtr_start'] = aln_start
    read_dtr_stats['dtr_end'] = aln_end

    return read_id, read_dtr_stats

def launch_dtr_pool( fastx, ftype, args ):
    
    func_args = []

    if ftype=="fastq":
        for read_num, (read_name, seq, qual) in enumerate(FastqGeneralIterator(open(fastx))):
            read_id = read_name.split(' ')[0]
            func_args.append( (read_id, seq) )

    elif ftype=="fasta":
        for read_num, (read_name, seq) in enumerate(SimpleFastaParser(open(fastx))):
            read_id = read_name.split(' ')[0]
            func_args.append( (read_id, seq) )
    
    results = launch_pool( args.threads, search_for_dtr, func_args )
    
    return dict(results)

def chunk_input_fastx(ftype, read_chunks, args):
    fastx_chunk_reads = defaultdict(list)
    for i,Seq in enumerate(SeqIO.parse(args.fastx, ftype)):
        chunk_id = i//args.chunk_size
        fastx_chunk_reads[chunk_id].append(Seq)
        
    fastx_chunk_fns = {}
    for i in range(len(read_chunks)):
        fastx_chunk_fns[i] = os.path.join(args.tmpdir, 'tmp.{}.{}'.format(i, ftype))
        SeqIO.write(fastx_chunk_reads[i], fastx_chunk_fns[i], ftype)

    return fastx_chunk_fns

def main(args):
    ftype  = check_input_format(args.fastx)
    n_reads = get_n_reads(args.fastx, ftype)
    os.makedirs(args.tmpdir, exist_ok=True)
    
    read_chunks = list(chunks(range(n_reads), args.chunk_size))

    fastx_chunk_fns = chunk_input_fastx(ftype, read_chunks, args)

    tmp_files = []
    for i,chunk in enumerate(tqdm(read_chunks)):

        dtr_stats = launch_dtr_pool( fastx_chunk_fns[i], \
                                     ftype,              \
                                     args)

        df_ = pd.DataFrame.from_dict(dtr_stats, orient='index').reset_index(drop=True)
        tmp_fname = os.path.join(args.tmpdir, 'chunk_{}.tmp'.format(i))
        tmp_files.append(tmp_fname)
        df_.to_csv(tmp_fname, sep='\t', index=False)

    df = pd.concat([pd.read_csv(fname, sep='\t') for fname in tmp_files], axis=0)

    if not args.no_hist:
        df[(df['has_dtr']==True) & (df['seq_len']<args.hist_max)].hist('seq_len', bins=50)
        plt.xlabel('DTR read length (bp)')
        plt.ylabel('Count')
        plt.savefig('{}.dtr.hist.png'.format(args.prefix))

    if not args.no_fasta:
        write_dtr_fasta(df, args, ftype)

    write_dtr_stats(df, args)

    shutil.rmtree(args.tmpdir)

if __name__=='__main__':
    args = parse_args()

    main(args)
