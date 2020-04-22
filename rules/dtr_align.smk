######################################################
# Align DTR seqs to ref to assess packaging strategy #
######################################################

rule align_dtr_seqs_to_polished_refs:
    input: 
        ref = BIN_CLUSTER_POLISHED_REF,
        binsdir = BINS_ROOT
    output: 
        plot = DTR_ALIGN_COORD_PLOT,
        tsv = DTR_ALIGN_TSV
    params: 
        prefix = DTR_ALIGN_PREFIX,
        ovlp = config['DTR_ALIGN']['ovlp']
    conda: '../envs/plot-umap.yml'
    threads: config['DTR_ALIGN']['threads']
    shell:
        'python {SCRIPT_DIR}/plot_headful_dtr_positions.py -t {threads} -p {params.prefix} '
        '-o {params.ovlp} {input.binsdir} {input.ref}'
