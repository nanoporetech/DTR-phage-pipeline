import os,sys
import shutil
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import argparse

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument('bins_dir', help='Fasta file containing read sequences', type=str)
    parser.add_argument('ref', help='Reference fasta that we align DTRs to', type=str)

    # Optional arguments
    parser.add_argument('-p', '--prefix', help='Output files prefix [output]', type=str, default='output')
    parser.add_argument('-o', '--overlap', help='Check for overlaps between the first and last <overlap> percent of the sequence [20]', type=int, default=20)
    parser.add_argument('-t', '--threads', help='Number of threads to use [4]', type=int, default=4)
    parser.add_argument('-c', '--chunk_size', help='Number of reads per chunk [5000]', type=int, default=5000)

    # Parse arguments
    args = parser.parse_args()

    # Get bin reads fasta based on the reference supplied
    args.bin_clust_id = args.ref.split('/')[-1].split('.')[0]
    args.bin_id = args.ref.split('/')[-1].split('_')[0]
    args.bin_fasta = os.path.join(args.bins_dir, args.bin_id, '{}.reads.fa'.format(args.bin_id))

    args.tmpdir = os.path.join('/'.join(args.prefix.split('/')[:-1]), 'dtr_tmp')

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

def call_minimap2(ref, query, paf):
    minimap2_CMD = 'minimap2 -t 1 -x map-ont --no-long-join -r100 {} {} > {} 2> /dev/null'.format(ref, query, paf)
    os.system(minimap2_CMD)

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
    call_minimap2(end_fn, start_fn, minimap2_paf)

    os.remove(start_fn)
    os.remove(end_fn)

    return minimap2_paf,ovlp_len

def get_n_reads(fasta):
    n_lines = 0
    with open(fasta) as f:
        for i, l in enumerate(f):
            n_lines += 1
    
    n_reads = len([read_tup for read_tup in SimpleFastaParser(open(fasta))])
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

    read_dtr_stats['seq_id'] = read_id
    read_dtr_stats['seq_len'] = len(seq)
    read_dtr_stats['has_dtr'] = dtr
    read_dtr_stats['dtr_length'] = aln_sum
    read_dtr_stats['dtr_start'] = aln_start
    read_dtr_stats['dtr_end'] = aln_end

    return read_id, read_dtr_stats

def fetch_dtr_seq(row, fasta_d):
    whole_seq = fasta_d[row['seq_id']]
    dtr_end = row['dtr_start'] + row['dtr_length']
    dtr_seq = str(whole_seq.seq[row['dtr_start']:dtr_end])
    record = SeqRecord(seq=Seq.Seq(dtr_seq), id='dtr_{}'.format(row['seq_id']), description='dtr_{}'.format(row['seq_id']))

    return dtr_seq, record

def launch_dtr_pool( fasta, args ):
    
    func_args = []

    for read_num, (read_name, seq) in enumerate(SimpleFastaParser(open(fasta))):
        read_id = read_name.split(' ')[0]
        func_args.append( (read_id, seq) )
    
    results = launch_pool( args.threads, search_for_dtr, func_args )
    
    return dict(results)

def chunk_input_fasta(read_chunks, args):
    fasta_chunk_reads = defaultdict(list)
    for i,Seq in enumerate(SeqIO.parse(args.bin_fasta, 'fasta')):
        chunk_id = i//args.chunk_size
        fasta_chunk_reads[chunk_id].append(Seq)
        
    fasta_chunk_fns = {}
    for i in range(len(read_chunks)):
        fasta_chunk_fns[i] = os.path.join(args.tmpdir, 'tmp.{}.fasta'.format(i))
        SeqIO.write(fasta_chunk_reads[i], fasta_chunk_fns[i], 'fasta')

    return fasta_chunk_fns

def align_dtr_seqs(df, args):
    fasta_d = SeqIO.index(args.bin_fasta, 'fasta')

    df['dtr_seq'] = ''
    seq_records = []
    for idx,row in df.iterrows():
        dtr_seq, seq_record = fetch_dtr_seq(row, fasta_d)
        seq_records.append(seq_record)

    dtr_seq_fname = os.path.join(args.tmpdir, 'dtr_seqs.fasta')
    SeqIO.write(seq_records, dtr_seq_fname, 'fasta')
    
    dtr_seq_paf = os.path.join(args.tmpdir, 'dtr_seqs.paf')
    call_minimap2(args.ref, dtr_seq_fname, dtr_seq_paf)
    
    dtr_df = pd.read_csv(dtr_seq_paf, sep='\t', usecols=range(12), names = ['qname', 'qlen', 'qstart', \
                                                                            'qend', 'strand', 'tname', \
                                                                            'tlen', 'tstart', 'tend',  \
                                                                            'matches', 'alnlen', 'mapqv'])
    return dtr_df

def main(args):
    n_reads = get_n_reads(args.bin_fasta)
    os.makedirs(args.tmpdir, exist_ok=True)
    
    read_chunks = list(chunks(range(n_reads), args.chunk_size))

    fasta_chunk_fns = chunk_input_fasta(read_chunks, args)

    tmp_files = []
    for i,chunk in enumerate(read_chunks):

        dtr_stats = launch_dtr_pool( fasta_chunk_fns[i], args)

        df_ = pd.DataFrame.from_dict(dtr_stats, orient='index').reset_index(drop=True)
        tmp_fname = os.path.join(args.tmpdir, 'chunk_{}.tmp'.format(i))
        tmp_files.append(tmp_fname)
        df_.to_csv(tmp_fname, sep='\t', index=False)

    df = pd.concat([pd.read_csv(fname, sep='\t') for fname in tmp_files], axis=0)

    dtr_df = align_dtr_seqs(df, args)

    dtr_df = dtr_df.sort_values('tend').reset_index(drop=False).rename(columns={'index':'orig_idx'}).reset_index()
    fig = plt.figure(figsize=[6,10])
    ax = fig.add_subplot(111)
    for idx,row in dtr_df.iterrows():
        ax.plot([row['tstart'], row['tend']], [idx, idx], linewidth=1)
    ax.set_xlabel('Position in {} reference (bp)'.format(args.bin_clust_id))
    ax.set_ylabel('All DTR seqs from bin {}'.format(args.bin_id))
    fig.savefig('{}.dtr.aligns.png'.format(args.prefix))

    dtr_df = dtr_df.loc[:,['qname','tstart','tend','tlen']]
    # Make column with tend as a fraction of total tlen
    dtr_df['tend_frac'] = dtr_df['tend'] / dtr_df['tlen']
    dtr_df['tend_frac'] = dtr_df['tend_frac'].round(4)
    dtr_df['bin_clust_id'] = args.bin_clust_id
    dtr_df.to_csv('{}.dtr.aligns.tsv'.format(args.prefix), sep='\t', index=False)

    shutil.rmtree(args.tmpdir)

if __name__=='__main__':
    args = parse_args()

    main(args)