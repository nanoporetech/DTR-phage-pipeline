#################################
# Create k-mer 2D embedding     #
#################################

rule calc_kmer_freq:
    input: FILTERED_FASTA
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
        fasta = FILTERED_FASTA
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
        fastx = FILTERED_FASTA
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
    """ Determines the bins """
    input: KMER_FREQS_UMAP
    output: KMER_BINS_MEMBERSHIP
    conda: '../envs/umap.yml'
    params:
        min_cluster = config['UMAP']['bin_min_reads'],
    shell:
        'python {SCRIPT_DIR}/run_hdbscan.py -o {output} -c {params.min_cluster} {input}'

rule list_kmer_bin_ids:
    """ 
    pull out list of bin IDs for quick parsing 
    we could checkpoint KMER_BINS_MEMBERSHIP instead, but this will be faster to parse each time we need it.
    """
    input: KMER_BINS_MEMBERSHIP
    output: KMER_BINS_LIST
    shell: "cut -f 6 {input} | tail -n+2 | sort | uniq > {output}"

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

############################################################################
# Generate dir for each k-mer bin                                          #       
#  - create read list                                                      #
#  - Get statistics for each kmer bin (genome size, coverage, etc)         #
############################################################################

checkpoint generate_bins_and_stats:
    """ 
    collects all bins stats into KMER_BIN_STATS,
    and, for each bin_id, generates in BINS_ROOT:
     * BIN_READLIST
     * BIN_GENOME_SIZE
    """
    input: 
        bin_membership=KMER_BINS_MEMBERSHIP,
    output: 
        all_stats=KMER_BIN_STATS,
        bins_dir=directory(BINS_ROOT)
    run:
        df_reads = pd.read_csv(input.bin_membership, sep='\t')

        # aggregate read stats into bin stats
        df_bins = df_reads.groupby('bin_id').agg(bases=pd.NamedAgg('length', sum),
                                                 genomesize=pd.NamedAgg('length', np.mean),
                                                 coverage=pd.NamedAgg('length', len),
                                                 readlist=pd.NamedAgg('read', set),
                                                )
        df_bins['rl_gs_est'] = True

        # write bin stats (without readlist) to stats file
        df_bins[['bases','genomesize','coverage','rl_gs_est']] \
               .sort_index() \
               .to_csv(output.all_stats, sep='\t')

        # write bin specific files
        for bin_id, row in df_bins.iterrows():
            #
            # create bin specific subdir of BINS_ROOT
            bin_dir = str(BIN_DIR).format(bin_id=bin_id)
            os.makedirs(bin_dir, exist_ok=True)
            #
            # save genome size to file
            gsize_fn = str(BIN_GENOME_SIZE).format(bin_id=bin_id)
            with open(gsize_fn, 'w') as fh:
                fh.write(f'{row["genomesize"]}\n')
            #
            # save read list to file
            rlist_fn = str(BIN_READLIST).format(bin_id=bin_id)
            with open(rlist_fn, 'w') as fh:
                fh.write('{}\n'.format('\n'.join(row['readlist'])))

def expand_template_from_bins(wildcards, template):
    # get dir through checkpoints to throw Exception if checkpoint is pending
    checkpoint_dir = checkpoints.generate_bins_and_stats.get(**wildcards).output
    # get bins from files
    bins, _bins = glob_wildcards(BIN_READLIST)
    # skips from config
    bins = [b for b in bins if b not in SKIP_BINS]
    # expand template
    return expand(str(template), bin_id=bins)

rule kmer_binned_fasta:
    input:
        read_list=BIN_READLIST,
        reads_fasta=FILTERED_FASTA,
    output: BIN_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.read_list} -o {output} {input.reads_fasta}'
