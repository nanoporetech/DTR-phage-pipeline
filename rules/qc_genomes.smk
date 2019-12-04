############################
# QC and summarize genomes #
############################

rule remove_adapters_from_ref_reads:
    input: ALL_POLISHED_COMBINED
    output: ALL_POLISHED_COMBINED_TRIMMED
    params:
        switches=config['PORECHOP']['switches']
    threads: config['PORECHOP']['threads']
    conda: '../envs/porechop.yml'
    shell:
        'porechop {params.switches} -t {threads} -i {input} '
        '-v 1 -o {output}'

rule write_dtr_report:
    input: ALL_POLISHED_COMBINED_TRIMMED
    output: ALL_POLISHED_COMBINED_TRIMMED_DTR_STATS
    params:
        polish_dir=BIN_CLUSTER_DIR
    conda: '../envs/mummer.yml'
    shell:
        'python {SCRIPT_DIR}/check_ref_for_DTR.py -o {output} {input} {params.polish_dir}'

rule build_bin_cluster_summary_table:
    input:
        dtr_table=ALL_POLISHED_COMBINED_TRIMMED_DTR_STATS,
    output: ALL_POLISHED_COMBINED_TRIMMED_DTR_GC_ABUND_STATS
    params:
        bins_dir=BINS_DIR,
        clust_dir=BIN_CLUSTER_DIR,
        prefix='{}{}'.format(SAMPLE,STYPE)
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/calculate_bin_cluster_abundance.py -p {params.prefix} '
        '-o {output} {input.dtr_table} {params.bins_dir} {params.clust_dir}'

rule combine_all_draft_stats:
    input:
        dtr_abund_stats=ALL_POLISHED_COMBINED_TRIMMED_DTR_GC_ABUND_STATS,
        cds_stats=ALL_POLISHED_COMBINED_CDS_SUMMARY,
        strand_stats=ALL_POLISHED_COMBINED_STRANDS
    output: ALL_POLISHED_STATS
    run:
        import pandas as pd
        df_gc     = pd.read_csv(input.dtr_abund_stats, sep='\t').rename(columns={'cluster_id':'cluster'}).set_index('cluster')
        df_cds    = pd.read_csv(input.cds_stats, sep='\t').rename(columns={'clust':'cluster'}).set_index('cluster')
        df_strand = pd.read_csv(input.strand_stats, sep='\t').rename(columns={'clust_id':'cluster','ref_read':'read'}).set_index('cluster')
        df_cds = df_cds.drop(['read','ref_len'], axis=1)
        df_strand['read'] = df_strand['read'].apply(lambda x: '_'.join(x.split('_')[:-2]))

        df_tmp = df_cds.join(df_strand, how='left')
        df     = df_tmp.join(df_gc, how='left').sort_values('cluster').reset_index()
        df.to_csv(output[0], sep='\t', index=False)

rule plot_all_prodigal_stats:
    input: ALL_POLISHED_STATS
    output:
        plot1=ALL_POLISHED_COMBINED_CDS_PLOT_ALL,
        plot2=ALL_POLISHED_COMBINED_CDS_PLOT_DTR,
        plot3=ALL_POLISHED_COMBINED_CDS_PLOT_DTR_NPOL10,
        plot4=ALL_POLISHED_COMBINED_CDS_PLOT_DTR_NPOL20,
    params:
        pol_dir=MEDAKA_DIR
    conda: '../envs/clustering.yml'
    shell:
        'python {SCRIPT_DIR}/combined_cds_summary_plot.py --output1={output.plot1} '
        '--output2={output.plot2} --output3={output.plot3} --output4={output.plot4} '
        '{input}'

rule create_dtr_fasta:
    input: 
        stats=ALL_POLISHED_COMBINED_TRIMMED_DTR_STATS,
        fasta=ALL_POLISHED_COMBINED_TRIMMED
    output: ALL_POLISHED_COMBINED_TRIMMED_DTR
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/write_dtr_refs_fasta.py -o {output} {input.stats} '
        '{input.fasta}'

################################################
# Check align strand of pol reads vs. ref read #
################################################

rule align_bin_cluster_reads_to_polished_ref:
    input: 
        ref=BIN_CLUSTER_POLISHED_REF,
        reads=BIN_CLUSTER_POL_READS_FASTA
    output: BIN_CLUSTER_POLISHED_POL_VS_REF_PAF
    threads: config['MINIMAP2']['STANDARD']['threads']
    params: 
        switches=config['MINIMAP2']['STANDARD']['switches']
    conda: '../envs/minimap2.yml'
    shell:
        'minimap2 -t {threads} {params.switches} {input.ref} '
        '{input.reads} > {output}'

rule parse_bin_cluster_reads_to_polished_ref_paf:
    input: 
        paf=BIN_CLUSTER_POLISHED_POL_VS_REF_PAF
    output: 
        counts=BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS,
        annots=BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS
    run:
        import os
        import re
        import pandas as pd
        clust_id = os.path.basename(input.paf).split('.')[0]
        cols = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', \
                'tlen', 'tstart', 'tend', 'matches', 'alnlen', 'mapqv']
        paf   = pd.read_csv(input.paf, sep='\t', header=None, usecols=list(range(12)), names=cols)
        paf   = paf.groupby(['qname','strand','tname']).size().reset_index(name='count')[['qname','strand','tname']]
        n_pos = paf[paf['strand']=='+'].shape[0]
        n_neg = paf[paf['strand']=='-'].shape[0]
        n_tot = n_pos + n_neg
        ref   = paf['tname'].unique()[0].split(':')[0]
        f_n   = open(output.counts, 'w')
        f_n.write('ref_read\tclust_id\tn_tot\tfrac_pos\tfrac_neg\n')
        f_n.write('{ref}\t{clust}\t{n_tot}\t{pos}\t{neg}\n'.format(ref=ref,clust=clust_id,n_tot=(n_pos+n_neg),pos=float(n_pos)/n_tot,neg=float(n_neg)/n_tot))
        
        reads_df       = paf[['qname','strand']].rename(columns={'qname':'read_id', 'strand':'label'})
        ref_df         = pd.DataFrame.from_dict({'read_id':[ref], 'label':['ref']})
        df             = pd.concat([reads_df, ref_df])
        df['clust_id'] = clust_id
        df.to_csv(output.annots, sep='\t', index=False)

rule combine_bin_cluster_strand_counts_into_table:
    input: 
        pol_dir=MEDAKA_DIR
    output:
        counts=ALL_POLISHED_COMBINED_STRANDS,
        annots=ALL_POLISHED_COMBINED_STRAND_ANNOTS
    run:
        import os
        import pandas as pd
        from glob import glob
        count_fns = glob(os.path.join(input.pol_dir, '*', '*.ref_read.strands.tsv'))
        count_dfs = [pd.read_csv(fn, sep='\t') for fn in count_fns]
        count_df  = pd.concat(count_dfs)
        count_df['bin']     = count_df['clust_id'].map(lambda x: int(x.split('_')[0]))
        count_df['cluster'] = count_df['clust_id'].map(lambda x: int(x.split('_')[1]))
        count_df            = count_df.sort_values(['bin','cluster']).drop(['bin','cluster'], axis=1).round({'frac_pos': 2, 'frac_neg': 2})
        count_df.to_csv(output.counts, sep='\t', index=False)
        annot_fns = glob(os.path.join(input.pol_dir, '*', '*.ref_read.strands.annots.tsv'))
        annot_dfs = [pd.read_csv(fn, sep='\t') for fn in annot_fns]
        annot_df  = pd.concat(annot_dfs)
        annot_df.to_csv(output.annots, sep='\t', index=False)

rule annotate_umap_with_strand_info:
    input:
        umap=KMER_FREQS_UMAP,
        annot=ALL_POLISHED_COMBINED_STRAND_ANNOTS
    output: ALL_POLISHED_COMBINED_STRAND_UMAP
    conda:'../envs/clustering.yml'
    shell:
        'python {SCRIPT_DIR}/plot_umap_with_strand_labels.py -o {output} '
        '{input.umap} {input.annot}'