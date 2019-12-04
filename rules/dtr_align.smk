######################################################
# Align DTR seqs to ref to assess packaging strategy #
######################################################

rule align_dtr_seqs_to_polished_refs:
    input: 
        ref = BIN_CLUSTER_POLISHED_REF
    output: 
        hist = DTR_ALIGN_HIST_PLOT,
        plot = DTR_ALIGN_COORD_PLOT,
        tsv = DTR_ALIGN_TSV
    params: 
        prefix = DTR_ALIGN_PREFIX,
        binsdir = BINS_DIR,
        ovlp = config['DTR_ALIGN']['ovlp']
    conda: '../envs/plot-umap.yml'
    threads: config['DTR_ALIGN']['threads']
    shell:
        'python {SCRIPT_DIR}/plot_headful_dtr_positions.py -t {threads} -p {params.prefix} '
        '-o {params.ovlp} {params.binsdir} {input.ref}'

rule combine_dtr_aligns:
    input:
        medaka_dir = MEDAKA_DIR
    output: 
        cyc_perm_stats = DTR_ALIGN_CYC_PERM_TSV
    run:
        import os
        import pandas as pd
        from glob import glob
        fns = glob(os.path.join(input.medaka_dir, '*', '*.dtr.aligns.tsv'))
        df = pd.concat([pd.read_csv(fn, sep='\t') for fn in fns])
        df['dtr_len'] = df['tend'] - df['tstart']
        df['left_dist'] = df['tend']
        df['right_dist'] = df['tlen'] - df['tstart']
        for_df = []
        for group,df_ in df.groupby('bin_clust_id'):
            cyc_perm = False
            # If the distance from the genome end to the inner DTR boundary is more than 2x the median DTR length, keep.
            # If > 25% of the reads pass that filter, label the cluster as being cyclically permuted.
            df_filt = df_[df_['left_dist'] > (2 * df_['dtr_len'].median())]
            df_filt = df_filt[df_filt['right_dist'] > (2 * df_filt['dtr_len'].median())]
            if df_filt.shape[0]>(0.25*df_.shape[0]):
                cyc_perm = True
            for_df.append( (df_.loc[0,'bin_clust_id'], cyc_perm))
        df_out = pd.DataFrame.from_records(for_df, columns=['cluster','cyc_perm'])
        df_out.to_csv(output.cyc_perm_stats, sep='\t', index=False)