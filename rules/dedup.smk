#######################
# Deduplicate genomes #
#######################

rule nucmer_for_duplicate_sequences:
    input: ALL_POLISHED_COMBINED_TRIMMED_DTR
    output: ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_DELTA
    params:
        prefix=ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_PREFIX,
    conda: '../envs/mummer.yml'
    shell:
        'nucmer -p {params.prefix} --nosimplify {input} {input}'

rule show_coords_for_duplicate_sequences:
    input: 
        delta=ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_DELTA,
    output:
        nucfilt=ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_COORDS
    conda: '../envs/mummer.yml'
    shell:
        'show-coords -c -l {input.delta} > {output.nucfilt}'

rule rewrite_dtr_fasta_without_seq_duplicates:
    input:
        fasta=ALL_POLISHED_COMBINED_TRIMMED_DTR,
        coords=ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_COORDS,
        stats=ALL_POLISHED_STATS,
    output: 
        fasta=ALL_POLISHED_COMBINED_TRIMMED_DTR_UNIQ,
        stats=ALL_POLISHED_STATS_FILT
    params:
        minpct=config['NUCMER']['minpct'],
        mincov=config['NUCMER']['mincov'],
        minpol=config['NUCMER']['minpol']
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/dedup_sample_genomes.py -p {params.minpct} '
        '-c {params.mincov} -n {params.minpol} -f {output.fasta} '
        '-s {output.stats} {input.coords} {input.fasta} {input.stats}'
