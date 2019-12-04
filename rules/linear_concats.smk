####################################
# Evaluate linear concatemer reads #
####################################

rule find_concatemer_read_ids:
    input: READS_IMPORT_FASTA
    output: CONCATEMER_READ_INFO
    params: 
        ovlp_len=config['LIN_CONCAT']['ovlp_len'],
        aln_dir=CONCATEMER_ALIGN_TMP_DIR,
        minlen=config['LIN_CONCAT']['minlen']
    threads: config['max_threads']
    conda: '../envs/plot-umap.yml'
    shell:
        'mkdir -p {params.aln_dir}; python {SCRIPT_DIR}/check_seqs_for_concatemers.py '
        '-o {params.ovlp_len} -d {params.aln_dir} -t {threads} -m {params.minlen} '
        '-o {output} {input}; '
        'rm -rf {params.aln_dir}'

rule plot_concatemer_copies_pdf:
    input: CONCATEMER_READ_INFO
    output: 
        copies=CONCATEMER_READ_PDF_PLOT_COPIES,
        lengths=CONCATEMER_READ_PDF_PLOT_LENGTHS,
    conda: '../envs/plot-umap.yml'
    shell:
        'python {SCRIPT_DIR}/plot_linear_concat_copy_numbers.py -c {output.copies} '
        '-l {output.lengths} {input}'

rule plot_concatemer_copy_length_contours:
    input: CONCATEMER_READ_INFO
    output: CONCATEMER_READ_COPY_REPEATS_CONTOURS
    params:
        sample=SAMPLE
    conda: '../envs/plot-umap.yml'
    shell:
        'python {SCRIPT_DIR}/plot_linear_concat_reads_hexbin.py -o {output} {input}'

rule grep_concatemer_seqs:
    input: 
        info=CONCATEMER_READ_INFO,
        fasta=READS_IMPORT_FASTA
    output: CONCATEMER_READ_FASTA
    conda: '../envs/seqkit.yml'
    shell:
        r"""grep -v read {input.info} | cut -d$'\t' -f1 | seqkit grep -f - {input.fasta} > {output}"""