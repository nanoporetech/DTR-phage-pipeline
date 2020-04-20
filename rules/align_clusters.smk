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

checkpoint combine_bin_cluster_read_info:
    input:
        clust_info_csvs=lambda w: expand(str(ALN_CLUST_OUTPUT_INFO), \
                                         bin_id=get_kmer_bins_good(w)),
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

def get_bin_clusters(wildcards):
    """ Use combo output from the checkpoint to get list of bin_clust_ids """
    chkpt = checkpoints.combine_bin_cluster_read_info
    with chkpt.get().output.combined.open() as table:
        return [f"{b}_{c}" \
                for b,c in pd.read_csv(table) \
                             .groupby(['bin_id', 'cluster']) \
                             .groups]

def get_polished_bin_outputs(wildcards, templates=POLISHED_TEMPLATE_LIST):
    bin_output_list = []
    final_bin_list = get_bin_clusters(wildcards)
    for filename_template in templates:
        bin_output_list.extend(expand(str(filename_template),
                                      bin_clust_id=final_bin_list))
    return bin_output_list

rule create_bin_cluster_read_lists:
    """ split bin cluster info by clusters """
    input: 
        # get bin_id form bin_clust_id
        clust_info=lambda w: str(ALN_CLUST_OUTPUT_INFO).format(
                                    bin_id=w.bin_clust_id.split("_")[0])
    output: 
        readlist=BIN_CLUSTER_READS_LIST,
        readinfo=BIN_CLUSTER_READS_INFO,
    run: 
        # parse table; filter by cluster
        bin_id, clust_id = wildcards.bin_clust_id.split("_")
        df_clust = pd.read_csv(input.clust_info).query('cluster == @clust_id')
        df_clust['read_id'].to_csv(output.readlist, index=False, header=False)
        df_clust.to_csv(output.readinfo, index=False)

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
