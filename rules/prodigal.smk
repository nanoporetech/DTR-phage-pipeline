######################################
# Prodigal on polished draft genomes #
######################################

rule prodigal:
    input: BIN_CLUSTER_POLISHED_REF 
    output: 
        fasta=BIN_CLUSTER_POLISHED_REF_PRODIGAL,
        txt=BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT
    conda: '../envs/prodigal-2.6.3.yml'
    threads: config['PRODIGAL']['threads'] 
    shell:
        'prodigal -p meta -q -i {input} -o {output.txt} -d {output.fasta}'

rule prodigal_summarize_result:
    input: 
        cds=BIN_CLUSTER_POLISHED_REF_PRODIGAL,
        ref=BIN_CLUSTER_POLISHED_REF
    output: BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/summarize_cds_result_fasta.py  -o {output} {input.cds} '
        '{input.ref}'

rule summarize_all_prodigal:
    output:
        table=ALL_POLISHED_COMBINED_CDS_SUMMARY,
    params:
        pol_dir=MEDAKA_DIR
    conda: '../envs/clustering.yml'
    shell:
        'python {SCRIPT_DIR}/combined_cds_summary.py -o {output.table} {params.pol_dir}'
