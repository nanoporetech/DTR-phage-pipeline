####################################################################################
# Kaiju Taxonomic Classification                                                   #
####################################################################################

rule kaiju:
    input: DTR_READS_FASTA
    output: KAIJU_RESULTS
    conda: '../envs/kaiju.yml'
    threads: config['max_threads']
    params:
        db = config['KAIJU']['db'],
        switch = config['KAIJU']['switch']
    shell:
        'kaiju {params.switch} -z {threads} -t {params.db}/nodes.dmp '
        '-f {params.db}/kaiju_db.fmi -i {input} -o {output}'

rule kaiju_taxa:
    input: KAIJU_RESULTS
    output: KAIJU_TAXIDS
    conda: '../envs/csvkit.yml'
    shell:
        'csvcut -t -c 3 {input} | sort -nk1 | uniq > {output}; '
        'echo -e 'taxID' | cat - {output} > /tmp/out && mv /tmp/out {output}'

rule kaiju_attach_taxa:
    input: 
        results = KAIJU_RESULTS,
    output: KAIJU_RESULTS_TAXA
    conda: '../envs/kaiju.yml'
    params:
        db = config['KAIJU']['db'],
    shell:
        'kaiju-addTaxonNames -p -t {params.db}/nodes.dmp -n {params.db}/names.dmp -i {input.results} -o {output}'

rule kaiju_to_krona:
    input: 
        results = KAIJU_RESULTS,
        dummy = KAIJU_RESULTS_TAXA
    output:
        krona = KAIJU_RESULTS_KRONA,
        html = KAIJU_RESULTS_KRONA_HTML
    conda: '../envs/kaiju.yml'
    params:
        db = config['KAIJU']['db'],
    shell:
        'kaiju2krona -u -t {params.db}/nodes.dmp -n {params.db}/names.dmp '
        '-i {input.results} -o {output.krona}; ktImportText -o {output.html} {output.krona}'