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

rule get_prodigal_result_statistics:
    input: 
        cds=BIN_CLUSTER_POLISHED_REF_PRODIGAL,
        ref=BIN_CLUSTER_POLISHED_REF
    output: BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS
    conda: '../envs/python.yml'
    shell:
        'python {SCRIPT_DIR}/summarize_cds_result_fasta.py  -o {output} {input.cds} '
        '{input.ref}'

####################################################
# Check strand abundance in each alignment cluster #
####################################################

rule align_pol_reads_to_polished_ref:
    input: 
        ref=BIN_CLUSTER_POLISHED_REF,
        reads=BIN_CLUSTER_POL_READS_FASTA
    output: BIN_CLUSTER_POLISHED_POL_VS_REF_PAF
    threads: config['MINIMAP2']['STANDARD']['threads']
    params: 
        switches=config['MINIMAP2']['STANDARD']['switches']
    conda: '../envs/minimap2.yml'
    shell:
        'minimap2 -t {threads} {params.switches} {input.ref} '
        '{input.reads} > {output}'

rule parse_align_clust_strands_paf:
    input: 
        paf=BIN_CLUSTER_POLISHED_POL_VS_REF_PAF
    output: 
        counts=BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS,
        annots=BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS
    run:
        clust_id = wildcards.bin_clust_id
        cols = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', \
                'tlen', 'tstart', 'tend', 'matches', 'alnlen', 'mapqv']
        paf   = pd.read_csv(input.paf, sep='\t', header=None, usecols=list(range(12)), names=cols)
        paf   = paf.groupby(['qname','strand','tname']).size().reset_index(name='count')[['qname','strand','tname']]
        n_pos = paf[paf['strand']=='+'].shape[0]
        n_neg = paf[paf['strand']=='-'].shape[0]
        n_tot = n_pos + n_neg
        ref   = paf['tname'].unique()[0].split(':')[0]
        f_n   = open(output.counts, 'w')
        f_n.write('read\tclust_id\tn_tot\tfrac_pos\tfrac_neg\n')
        f_n.write('{ref}\t{clust}\t{n_tot}\t{pos}\t{neg}\n'.format(ref=ref,clust=clust_id,n_tot=(n_pos+n_neg),pos=float(n_pos)/n_tot,neg=float(n_neg)/n_tot))
        
        reads_df       = paf[['qname','strand']].rename(columns={'qname':'read_id', 'strand':'label'})
        ref_df         = pd.DataFrame.from_dict({'read_id':[ref], 'label':['ref']})
        df             = pd.concat([reads_df, ref_df])
        df['clust_id'] = clust_id
        df.to_csv(output.annots, sep='\t', index=False)
