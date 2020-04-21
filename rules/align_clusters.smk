#####################################
# All-vs-all alignment for each bin #
#####################################

rule bin_ava_minimap2:
    input: BIN_FASTA
    output: ALN_CLUST_PAF
    conda: '../envs/minimap2.yml'
    params:
        cli = config['BIN_AVA_ALIGN']['MINIMAP2']['cli']
    threads: config['BIN_AVA_ALIGN']['MINIMAP2']['threads']
    shell:
        'minimap2 -t {threads} {params.cli} {input} {input} > {output}'

rule bin_ava_clustering:
    input: ALN_CLUST_PAF,
    output:
        heatmap = ALN_CLUST_OUTPUT_HEATMAP,
        info = ALN_CLUST_OUTPUT_INFO
    conda: '../envs/clustering.yml'
    params:
        prefix = str(ALN_CLUST_OUTPUT_PREFIX),
        min_score_frac = config['BIN_AVA_ALIGN']['CLUST']['min_score_frac'],
        min_reads = config['BIN_AVA_ALIGN']['CLUST']['min_reads'],
    shell:
        'python {SCRIPT_DIR}/cluster_ava_alignments.py -p {params.prefix} '
        ' -s {params.min_score_frac} -n {params.min_reads} {input}'

rule combine_bin_cluster_read_info:
    input:
        clust_info_csvs=lambda w: expand_template_from_bins(w, \
                                                            ALN_CLUST_OUTPUT_INFO),
    output: 
        combined=ALN_CLUST_READS_COMBO
    run:
        fns = input.clust_info_csvs
        df = pd.concat([pd.read_csv(fn) for fn in fns], ignore_index=True)
        df = df.sort_values(['bin_id', 'cluster'])
        print(f"Writing to {output.combined}")
        df.to_csv(str(output.combined), index=False)

############################################
# Separate out read ID's from each cluster #
############################################

checkpoint create_bin_cluster_read_lists:
    """ split all bins into clusters """
    input: 
        clust_info=lambda w: expand_template_from_bins(w, ALN_CLUST_OUTPUT_INFO)
    output: 
        directory(BIN_CLUSTER_ROOT)
    run: 
        # parse table; filter by cluster
        for clust_info_file in input.clust_info:
            bin_id = os.path.basename(os.path.dirname(clust_info_file))
            print(f"DEBUG: bin id: {bin_id}")
            df = pd.read_csv(clust_info_file)
            for clust_id, df_clust in df.groupby('cluster'):
                bin_clust_id = f"{bin_id}_{clust_id}"
                bin_clust_dir = str(BIN_CLUSTER_DIR).format(bin_clust_id=bin_clust_id)
                os.makedirs(bin_clust_dir, exist_ok=True)
                readlist = str(BIN_CLUSTER_READS_LIST).format(bin_clust_id=bin_clust_id)
                df_clust['read_id'].to_csv(readlist, index=False, header=False)
                readinfo = str(BIN_CLUSTER_READS_INFO).format(bin_clust_id=bin_clust_id)
                df_clust.to_csv(readinfo, index=False)

def expand_template_from_bin_clusters(wildcards, template):
    """ looks for "{bin_id}_{vlust_id}" folders in BIN_CLUSTER_ROOT """
    # get dir through checkpoints to throw Exception if checkpoint is pending
    checkpoint_dir = checkpoints.create_bin_cluster_read_lists.get(**wildcards).output[0]
    # get bin_bin_clusters from files 
    #bin_clusters, = glob_wildcards(BIN_CLUSTER_READS_LIST)
    # use chkpt dir instead of template because somethings not working
    bin_clusters, = glob_wildcards(os.path.join(checkpoint_dir, \
                                                "{bin_clust_id}", \
                                                "readlist.csv"))

    print(f"DEBUG: clusters: {repr(bin_clusters)}")
    # expand template
    return expand(str(template), bin_clust_id=bin_clusters)

rule get_cluster_reads_fasta:
    input: 
        clustreads=BIN_CLUSTER_READS_LIST,
        allfasta=FILTERED_FASTA
    output: BIN_CLUSTER_READS_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.clustreads} {input.allfasta} > {output}'

rule select_read_with_longest_DTR:
    input: 
        cluster_info=BIN_CLUSTER_READS_INFO
    output: 
        ref=BIN_CLUSTER_REF_READ_LIST,
        pol=BIN_CLUSTER_POL_READS_LIST
    run:
        clust_info_df = pd.read_csv(input.cluster_info)
        ref_idx  = clust_info_df['clust_read_score'].idxmax()
        ref_read = clust_info_df.loc[ref_idx, ['read_id']]
        pol_reads = clust_info_df[~clust_info_df['read_id'].isin(ref_read.values)]['read_id']
        ref_read.to_csv(output.ref, index=False, header=False)
        pol_reads.to_csv(output.pol, index=False, header=False)

rule create_bin_cluster_ref_fasta:
    input: 
        readlist=BIN_CLUSTER_REF_READ_LIST,
        fasta=FILTERED_FASTA,
    output: BIN_CLUSTER_REF_READ_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.readlist} {input.fasta} > {output}'

rule create_bin_cluster_pol_fastqs:
    input: 
        readlist=BIN_CLUSTER_POL_READS_LIST,
        fasta=FILTERED_FASTA,
    output: BIN_CLUSTER_POL_READS_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.readlist} {input.fasta} > {output}'
