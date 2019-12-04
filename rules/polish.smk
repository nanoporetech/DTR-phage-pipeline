######################################################################
# Run multiple rounds of Racon on the ref read from each bin cluster #
######################################################################

rule racon_polishing:
    input: 
        # raw_reads=BIN_CLUSTER_POL_READS_FASTQ,
        raw_reads=BIN_CLUSTER_POL_READS_FASTA,
        draft=BIN_CLUSTER_REF_READ_FASTA
    output: 
        polished=BIN_CLUSTER_RACON_POLISHED_FASTA
    params:
        iterations=RACON_ROUNDS,
        minq=config['RACON']['quality-threshold']
    conda: '../envs/racon.yml'
    threads: config['RACON']['threads']
    shell:
        'python {SCRIPT_DIR}/run_racon.py -t {threads} -n {params.iterations} '
        '-q {params.minq} -o {output.polished} {input.raw_reads} {input.draft}'

#################################################
# Run Medaka on each bin cluster ref read       #
#################################################

rule run_medaka:
    input:
        raw_reads=BIN_CLUSTER_POL_READS_FASTA,
        draft=BIN_CLUSTER_RACON_POLISHED_FASTA
    output:
        BIN_CLUSTER_POLISHED_REF_TMP
    params:
        out_dir=lambda x: str(Path(BIN_CLUSTER_POLISHED_REF_TMP.format(**x)).parent),
        model=config['MEDAKA']['model']
    threads: config['MEDAKA']['threads']
    conda: '../envs/medaka-0.11.0.yml'
    shell:
        'medaka_consensus -i {input.raw_reads} -d {input.draft} '
        '-o {params.out_dir} -t {threads} -m {params.model}; '
        'mv {params.out_dir}/consensus.fasta {output}'

rule rename_polished_ref_reads:
    input: BIN_CLUSTER_POLISHED_REF_TMP
    output: BIN_CLUSTER_POLISHED_REF
    conda: '../envs/seqkit.yml'
    shell: 
         "seqkit replace -p ':\d+\.\d-\d+\.\d' -r '_{wildcards.bin_clust_id}' "
         '{input} > {output} '

rule combine_all_polished_ref_reads:
    input: 
        dirname=MEDAKA_DIR,
    output: ALL_POLISHED_COMBINED
    shell:
        'cat {input.dirname}/*/*.medaka.fasta > {output}'