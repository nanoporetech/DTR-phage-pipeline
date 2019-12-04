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
        prefix = ALN_CLUST_OUTPUT_PREFIX,
        min_score_frac = config['BIN_AVA_ALIGN']['CLUST']['min_score_frac'],
        min_reads = config['BIN_AVA_ALIGN']['CLUST']['min_reads'],
    shell:
        'python {SCRIPT_DIR}/find_full_len_align_clusts.py -p {params.prefix} '
        ' -s {params.min_score_frac} -n {params.min_reads} {input}'

rule combine_bin_cluster_read_info:
    input:
        read_select_dir=ALN_CLUST_DIR
    output: 
        combined=ALN_CLUST_READS_COMBO
    run:
        import os,sys
        import pandas as pd
        from glob import glob
        read_select_dir = str(input.read_select_dir)
        fns = glob(os.path.join(read_select_dir, '*', '*.clust.info.csv'))
        df = pd.concat([pd.read_csv(fn) for fn in fns], ignore_index=True)
        df = df.sort_values(['bin_id', 'cluster'])
        df.to_csv(output.combined, index=False)

############################################
# Separate out read ID's from each cluster #
############################################

rule create_bin_cluster_read_lists:
    input: 
        clust_reads=ALN_CLUST_READS_COMBO
    output:
        readlist=dynamic(BIN_CLUSTER_READS_LIST),
        readinfo=dynamic(BIN_CLUSTER_READS_INFO),
    params:
        bin_cluster_dir=BIN_CLUSTER_DIR
    run: 
        import os,sys
        import pandas as pd
        df   = pd.read_csv(input.clust_reads)
        df_g = df.groupby(['bin_id', 'cluster'])
        for (bin_id,clust),df_ in df_g:
            bin_clust_id = '{}_{}'.format(str(bin_id), str(clust))
            clust_outdir = Path(os.path.join(params.bin_cluster_dir, bin_clust_id))
            clust_outdir.mkdir(exist_ok=True)
            clust_reads_file = clust_outdir / '{}.readlist.csv'.format(bin_clust_id)
            df_['read_id'].to_csv(clust_reads_file, index=False)
            clust_info_file = clust_outdir / '{}.readinfo.csv'.format(bin_clust_id)
            df_.to_csv(clust_info_file, index=False)

rule get_cluster_reads_fasta:
    input: 
        clustreads=BIN_CLUSTER_READS_LIST,
        allfasta=DTR_READS_FASTA
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
        import pandas as pd
        clust_info_df = pd.read_csv(input.cluster_info)
        ref_idx  = clust_info_df['clust_read_score'].idxmax()
        ref_read = clust_info_df.loc[ref_idx, ['read_id']]
        pol_reads = clust_info_df[~clust_info_df['read_id'].isin(ref_read.values)]['read_id']
        ref_read.to_csv(output.ref, index=False)
        pol_reads.to_csv(output.pol, index=False)

rule create_bin_cluster_ref_fasta:
    input: 
        readlist=BIN_CLUSTER_REF_READ_LIST,
        fasta=DTR_READS_FASTA,
    output: BIN_CLUSTER_REF_READ_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.readlist} {input.fasta} > {output}'

rule create_bin_cluster_pol_fastqs:
    input: 
        readlist=BIN_CLUSTER_POL_READS_LIST,
        fasta=DTR_READS_FASTA,
    output: BIN_CLUSTER_POL_READS_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        'seqkit grep -f {input.readlist} {input.fasta} > {output}'