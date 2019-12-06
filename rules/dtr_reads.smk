####################################
# Characterize DTR reads in sample #
####################################

rule find_all_DTR_reads:
    input: READS_IMPORT_FASTA
    output:
        fasta = DTR_READS_FASTA,
        stats = DTR_READS_STATS,
        hist = DTR_READS_HIST,
    conda: '../envs/dtr.yml'
    params:
        prefix = DTR_READS_PREFIX,
        tmpdir = DTR_READS_TMP_DIR,
        ovlp   = config['DTR_FIND']['ovlp'],
        chunksize = config['DTR_FIND']['chunksize']
    threads: config['max_threads']
    shell:
        'python {SCRIPT_DIR}/find_dtr_all_seqs.py -t {threads} -p {params.prefix} '
        '-d {params.tmpdir} -o {params.ovlp} -c {params.chunksize} {input}'