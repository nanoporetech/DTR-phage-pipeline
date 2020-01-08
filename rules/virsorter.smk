####################################
# Identify phage reads in sample   #
####################################

# location of Virsorter output files
phage_fasta_files = expand( \
    "{dir}/Predicted_viral_sequences/VIRSorter_cat-{cat}.fasta", \
    dir=VIRSORTER_DIR, \
    cat=[1,2,3])
prophage_fasta_files = expand( \
    "{dir}/Predicted_viral_sequences/VIRSorter_vprophages_cat-{cat}.fasta", \
    dir=VIRSORTER_DIR, \
    cat=[4,5,6])

# which sequences to return
if config['VIRSORTER']['prophage']:
    # return phage and prophage
    fasta_list = phage_fasta_files + prophage_fasta_files
else:
    # return just the phage
    fasta_list = phage_fasta_files

rule virsorter_identify_phage:
    """ Run VIRSORTER """
    input: READS_IMPORT_FASTA
    output:
        table=VIRSORTER_DIR / "VIRSorter_global-phage-signal.csv",
        fasta=phage_fasta_files + prophage_fasta_files
    conda: '../envs/virsorter.yml'
    params:
        output_dir=VIRSORTER_DIR,
        db=config['VIRSORTER']['db'],
        virome="--virome" if config['VIRSORTER']['virome'] else "",
        diamond="--diamond" if config['VIRSORTER']['diamond'] else "",
    threads: config['max_threads']
    shell:
        'wrapper_phage_contigs_sorter_iPlant.pl -f {input} \
         --wdir {params.output_dir} \
         --data-dir {params.data_dir} \
         --db {params.db} --ncpu {threads} \
         {params.virome} {params.diamond}'

rule virsorter_collect_phage:
    """ Collect VIRSORTER output """
    input: fasta_list
    output:
        fasta=VIRSORTER_FASTA
    shell: 'cat {input} > {output}'
