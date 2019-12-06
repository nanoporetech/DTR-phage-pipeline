#######################
# Deduplicate genomes #
#######################

rule nucmer_for_duplicate_sequences:
    input: ALL_POL
    output: temp(ALL_POL_NUCMER_DELTA)
    params:
        prefix=ALL_POL_NUCMER_PREFIX,
    conda: '../envs/mummer.yml'
    shell:
        'nucmer -p {params.prefix} --nosimplify {input} {input}'

rule show_coords_for_duplicate_sequences:
    input: 
        delta=ALL_POL_NUCMER_DELTA,
    output:
        nucfilt=temp(ALL_POL_NUCMER_COORDS)
    conda: '../envs/mummer.yml'
    shell:
        'show-coords -c -l -T {input.delta} > {output.nucfilt}'

rule rewrite_pol_fasta_without_seq_duplicates:
    input:
        fasta=ALL_POL,
        coords=ALL_POL_NUCMER_COORDS,
        stats=ALL_POL_STATS,
    output: 
        fasta=ALL_POL_UNIQ,
        stats=ALL_POL_STATS_FILT
    params:
        minpct=config['NUCMER']['minpct'],
        mincov=config['NUCMER']['mincov'],
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/dedup_sample_genomes.py -p {params.minpct} '
        '-c {params.mincov} -f {output.fasta} -s {output.stats} '
        '{input.coords} {input.fasta} {input.stats}'
