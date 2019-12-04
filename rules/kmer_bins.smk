#################################
# Create k-mer 2D embedding     #
#################################

rule calc_kmer_freq:
    input: DTR_READS_FASTA
    output: temp(KMER_FREQS_TMP)
    conda: '../envs/python.yml'
    params: 
        k = config['KMER_CALC']['k'],
        cli = config['KMER_CALC']['cli']
    threads: config['max_threads']
    shell:
        'python {SCRIPT_DIR}/kmer_freq.py {params.cli} -k {params.k} -t {threads} '
        '{input} > {output}'

rule add_qscore_to_kmer_freq:
    input: 
        freq=KMER_FREQS_TMP,
        summ=READS_IMPORT_SUMMARY
    output: KMER_FREQS
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/add_qscore_to_kmer_freq.py -o {output} {input.freq} {input.summ}'

rule create_umap_kmer_freq_map:
    input: 
        comp = KMER_FREQS,
    output: KMER_FREQS_UMAP
    conda: '../envs/umap.yml'
    params:
        min_rl = config['UMAP']['min_rl'],
        min_q = config['UMAP']['min_q'],
        min_d =config['UMAP']['min_d'],
        neighbors=config['UMAP']['neighbors']
    shell:
        'python {SCRIPT_DIR}/run_umap.py -l {params.min_rl} -q {params.min_q} '
        '-d {params.min_d} -n {params.neighbors} -o {output} {input.comp}'

#################################
# Create k-mer map plots        #
#################################

rule label_umap_freq_map_with_qscore:
    input: KMER_FREQS_UMAP,
    output: KMER_FREQS_UMAP_QSCORE
    params:
        size = config['UMAP']['scatter_size'],
        alpha = config['UMAP']['scatter_alpha'],
    conda: '../envs/plot-umap.yml'
    shell:
        'python {SCRIPT_DIR}/plot_umap_with_qscore.py -o {output} -s {params.size} '
        '-a {params.alpha} {input}'

rule label_umap_freq_map_with_gc_content:
    input: 
        umap = KMER_FREQS_UMAP,
        fasta = DTR_READS_FASTA
    output: KMER_FREQS_UMAP_GC
    params:
        size = config['UMAP']['scatter_size'],
        alpha = config['UMAP']['scatter_alpha'],
    conda: '../envs/plot-umap.yml'
    shell:
        'python {SCRIPT_DIR}/plot_umap_with_gc_content.py -o {output} -s {params.size} '
        '-a {params.alpha} {input.umap} {input.fasta}'

rule label_umap_freq_map_with_readlength:
    input: 
        umap = KMER_FREQS_UMAP,
        fastx = DTR_READS_FASTA
    output: KMER_FREQS_UMAP_READLENGTH
    params:
        min_rl=config['UMAP']['min_rl'],
        max_cb_rl=config['UMAP']['max_cb_rl'],
        size = config['UMAP']['scatter_size'],
        alpha = config['UMAP']['scatter_alpha']
    conda: '../envs/plot-umap.yml'
    shell:
        'python {SCRIPT_DIR}/plot_umap_with_readlength.py -o {output} -m {params.min_rl} '
        '-n {params.max_cb_rl} -s {params.size} -a {params.alpha} {input.umap} {input.fastx}'

rule call_umap_freq_map_bins:
    input: KMER_FREQS_UMAP
    output: KMER_BINS_MEMBERSHIP
    conda: '../envs/umap.yml'
    params:
        min_cluster = config['UMAP']['bin_min_reads'],
    shell:
        'python {SCRIPT_DIR}/run_hdbscan.py -o {output} -c {params.min_cluster} {input}'

rule plot_umap_freq_map_bins:
    input: KMER_BINS_MEMBERSHIP
    output: KMER_FREQS_UMAP_BINS_PLOT
    conda: '../envs/plot-umap.yml'
    params:
        size = config['UMAP']['scatter_size'],
        alpha = config['UMAP']['scatter_alpha'],
    shell:
        'python {SCRIPT_DIR}/plot_umap_with_hdbscan_labels.py -o {output} -s {params.size} '
        '-a {params.alpha} {input}'

rule label_umap_freq_map_with_tax:
    input: 
        umap = KMER_FREQS_UMAP,
        annot = KAIJU_RESULTS_TAXA,
    output: KMER_FREQS_UMAP_TAX
    conda: '../envs/plot-umap.yml'
    params:
        size = config['UMAP']['scatter_size'],
        alpha = config['UMAP']['scatter_alpha'],
    shell:
        'python {SCRIPT_DIR}/plot_umap_with_kaiju_labels.py -r {wildcards.rank} -o {output} '
        '-s {params.size} -a {params.alpha} {input.umap} {input.annot}'

###################################
# Generate dir for each k-mer bin #
###################################

rule create_kmer_bins:
    input: 
        bin_membership = KMER_BINS_MEMBERSHIP,
        all_reads = KMER_FREQS,
    output: dynamic(BIN_READLIST)
    params:
        binned_dir = str(BINS_DIR)
    run:
        import os
        from pathlib import Path
        import pandas as pd
        df = pd.read_csv(input.bin_membership, sep='\t')
        for j in df['bin_id'].unique():
            outdir = Path(os.path.join(params.binned_dir, str(j)))
            outdir.mkdir(exist_ok=True)
            read_file = outdir / 'read_list.txt'
            with read_file.open('w') as fh:
                fh.write('%s\n' % '\n'.join(df[df['bin_id']==j]['read'].unique()))

rule kmer_binned_fasta:
    input:
        read_list=BIN_READLIST,
        reads_fasta=DTR_READS_FASTA,
    output: BIN_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.read_list} -o {output} {input.reads_fasta}'

####################################################################################
# Get statistics for each kmer bin (genome size, coverage, etc)                    #
####################################################################################

rule generate_bin_stats:
    input: 
        bins = KMER_BINS_MEMBERSHIP
    output: 
        stats = KMER_BIN_STATS
    params:
        bins_dir = str(BINS_DIR)
    run:
        import pandas as pd
        import os
        df = pd.read_csv(input.bins, sep='\t')
        df['bases'] = df.groupby('bin_id')['length'].transform('sum')
        df['genomesize'] = df.groupby('bin_id')['length'].transform('mean')
        df['coverage'] = df['bases'] / df['genomesize']
        df['rl_gs_est'] = True
        df = df.loc[:,['bin_id','bases','genomesize','coverage','rl_gs_est']].sort_values('bin_id').drop_duplicates()
        df.to_csv(output.stats, sep='\t', index=False)
        for idx,row in df.iterrows():
            fn = os.path.join(params.bins_dir, str(row['bin_id']), 'genome_size.txt')
            with open(fn, 'w') as f:
                f.write('{}\n'.format(row['genomesize']))