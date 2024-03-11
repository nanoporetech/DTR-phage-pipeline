from pathlib import Path
import os
from glob import glob
import pandas as pd
import numpy as np
from snakemake.logging import logger

configfile: 'config.yml'
# this SCRIPT_DIR definition uses something outside the official API
# snakemake may change out from under us.
# To reduce the chance of future problems, we should switch all the scripts
# to snakemake scripts from shell scripts....
SCRIPT_DIR = workflow.basedir + "/scripts"    

###################################
# INPUT VARS                      #
###################################
SAMPLE        = config['sample'] 
STYPE         = config['stype']
K             = config['KMER_CALC']['k']
DATABASE_NAME = ["nr_euk"] # can also point to dbs like "progenomes" and "mar"
TAX_RANK      = ["0","1","2","3","4","5","6"]


###################################
# OUTPUT VARS                     #
###################################
OUTPUT_ROOT  = config['output_root']
if OUTPUT_ROOT.endswith('/'):
    OUTPUT_ROOT = OUTPUT_ROOT[:-1]
OUTPUT_ROOT  = Path(OUTPUT_ROOT)
VERSION      = config['version']
RACON_ROUNDS = config['RACON']['repeat']


###################################
# READ IMPORT FILES               #
###################################
READS_IMPORT_FASTA      = config['input_fasta']
READS_IMPORT_SUMMARY    = config['input_summary']


###################################
# SUMMARY STATS FILES             #
###################################
READS_DIR      = OUTPUT_ROOT / SAMPLE / STYPE / VERSION / 'reads_summary'
SUMMARY_PLOT   = READS_DIR / 'reads.summary.stats.png'


###################################
# DTR FINDING FILES               #
###################################
DTR_READS_DIR                = OUTPUT_ROOT / SAMPLE / STYPE / VERSION / 'dtr_reads'
DTR_READS_PREFIX             = '{}/output'.format(str(DTR_READS_DIR))
DTR_READS_TMP_DIR            = DTR_READS_DIR / 'aln_tmp'
DTR_READS_FASTA              = '{}.dtr.fasta'.format(str(DTR_READS_PREFIX))
DTR_READS_STATS              = '{}.dtr.stats.tsv'.format(str(DTR_READS_PREFIX))
DTR_READS_HIST               = '{}.dtr.hist.png'.format(str(DTR_READS_PREFIX))

###################################
# VIRSORTER FILES                 #
###################################
VIRSORTER_DIR            = OUTPUT_ROOT / SAMPLE / STYPE / VERSION / 'virsorter'
VIRSORTER_FASTA          = VIRSORTER_DIR / 'all_phage.fasta'

pre_filter = config.get('pre_filter', 'DTR').upper()
FILTERED_FASTA = VIRSORTER_FASTA if pre_filter == "VIRSORTER" \
                                 else DTR_READS_FASTA if pre_filter == "DTR" \
                                 else READS_IMPORT_FASTA


###################################
# KAIJU CLASSIFICATION            #
###################################
KAIJU_DB_DIR              = config['KAIJU']['db']
KAIJU_DIR                 = OUTPUT_ROOT /  SAMPLE / STYPE / VERSION / 'kaiju'
KAIJU_RESULTS             = KAIJU_DIR / 'results.tsv'
KAIJU_RESULTS_TAXA        = KAIJU_DIR / 'results.taxa.tsv'
KAIJU_RESULTS_KRONA       = KAIJU_DIR / 'results.krona'
KAIJU_RESULTS_KRONA_HTML  = KAIJU_DIR / 'results.html'
KAIJU_TAXIDS              = KAIJU_DIR / 'taxids.csv'
KAIJU_TAXINFO             = KAIJU_DIR / 'taxinfo.csv'
KAIJU_TAX_GENOMESIZE_EST  = KAIJU_DIR / 'taxid_genomesize_estimates.csv'


###################################
# K-mer UAMP files and plots      #
###################################
KMER_BIN_ROOT                     = OUTPUT_ROOT /  SAMPLE / STYPE / VERSION / 'kmer_binning'
KMER_BINS_MEMBERSHIP              = KMER_BIN_ROOT / 'bin_membership.tsv'
KMER_BINS_LIST                    = KMER_BIN_ROOT / 'bin_list.txt'
KMER_BIN_STATS                    = KMER_BIN_ROOT / 'bin_stats.csv'
KMER_FREQS_TMP                    = KMER_BIN_ROOT / 'kmer_comp.tmp'
KMER_FREQS                        = KMER_BIN_ROOT / 'kmer_comp.tsv'
KMER_FREQS_UMAP                   = KMER_BIN_ROOT / 'kmer_comp.umap.tsv'
KMER_FREQS_UMAP_TAX               = KMER_BIN_ROOT / 'kmer_comp.umap.{database}.{rank}.png'
KMER_FREQS_UMAP_QSCORE            = KMER_BIN_ROOT / 'kmer_comp.umap.qscore.png'
KMER_FREQS_UMAP_GC                = KMER_BIN_ROOT / 'kmer_comp.umap.gc.png'
KMER_FREQS_UMAP_READLENGTH        = KMER_BIN_ROOT / 'kmer_comp.umap.readlength.png'
KMER_FREQS_UMAP_BINS_PLOT         = KMER_BIN_ROOT / 'kmer_comp.umap.bins.png'
KMER_FREQS_UMAP_BINS_COORDS       = KMER_BIN_ROOT / 'kmer_comp.umap.bins.tsv'


###################################
# BIN ANALYSIS FILES              #
###################################
BINS_ROOT                   = KMER_BIN_ROOT / 'bins'
BIN_DIR                     = BINS_ROOT / '{bin_id}'
BIN_READLIST                = BIN_DIR / 'read_list.txt'
BIN_FASTA                   = BIN_DIR / '{bin_id}.reads.fa'
BIN_GENOME_SIZE             = BIN_DIR / 'genome_size.txt'
SKIP_BINS                   = set(str(b) for b in \
                                  config['SKIP_BINS'][SAMPLE][STYPE][VERSION])
BINNED_ANALYSIS_ROOT        = KMER_BIN_ROOT / "refine_bins"

#####################################
# Alignment clustering  FILES       #
#####################################
ALN_CLUST_DIR            = BINNED_ANALYSIS_ROOT / 'alignments'
ALN_CLUST_PAF            = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.ava.paf'
ALN_CLUST_OUTPUT_PREFIX  = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust'
ALN_CLUST_OUTPUT_HEATMAP = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust.heatmap.png'
ALN_CLUST_OUTPUT_INFO    = ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust.info.csv'
ALN_CLUST_READS_COMBO    = ALN_CLUST_DIR / 'all_bins.clust.info.csv'


###################################
# Separate cluster-specific reads #
###################################
BIN_CLUSTER_ROOT                = BINNED_ANALYSIS_ROOT / 'align_cluster_reads'
BIN_CLUSTER_DIR                 = BIN_CLUSTER_ROOT / '{bin_clust_id}'
BIN_CLUSTER_READS_INFO          = BIN_CLUSTER_DIR / '{bin_clust_id}.readinfo.csv'
BIN_CLUSTER_READS_LIST          = BIN_CLUSTER_DIR / 'readlist.csv'
BIN_CLUSTER_READS_FASTA         = BIN_CLUSTER_DIR / '{bin_clust_id}.reads.fasta'
BIN_CLUSTER_REF_READ_LIST       = BIN_CLUSTER_DIR / '{bin_clust_id}.ref_readlist.csv'
BIN_CLUSTER_POL_READS_LIST      = BIN_CLUSTER_DIR / '{bin_clust_id}.pol_readlist.csv'
BIN_CLUSTER_REF_READ_FASTA      = BIN_CLUSTER_DIR / '{bin_clust_id}.ref_read.fa'
BIN_CLUSTER_POL_READS_FASTQ     = BIN_CLUSTER_DIR / '{bin_clust_id}.pol_reads.fq'
BIN_CLUSTER_POL_READS_FASTA     = BIN_CLUSTER_DIR / '{bin_clust_id}.pol_reads.fa'


####################################
# Polish cluster reads using Racon #
####################################
POLISH_DIR                       = BINNED_ANALYSIS_ROOT / 'align_cluster_polishing'
RACON_DIR                        = POLISH_DIR / 'racon'
BIN_CLUSTER_RACON_POLISHED_FASTA = RACON_DIR / '{bin_clust_id}' / '{{bin_clust_id}}.ref_read.racon_{repeats}x.fasta'.format(repeats=RACON_ROUNDS)


#####################################
# Polish cluster reads using Medaka #
#####################################
MEDAKA_DIR                                    = POLISH_DIR / 'medaka'
BIN_CLUSTER_POLISHED_REF_TMP                  = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.tmp.fasta'
BIN_CLUSTER_POLISHED_REF                      = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.fasta'
BIN_CLUSTER_POLISHED_POL_VS_REF_PAF           = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.paf'
BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS       = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.strands.summary.tsv'
BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.strands.reads.tsv'
BIN_CLUSTER_POLISHED_REF_PRODIGAL             = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.fasta'
BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT         = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.txt'
BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS       = MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.stats.txt'

#########################################
# Check for fixed or cyclic permutation #
#########################################
DTR_ALIGN_PREFIX                 = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read')
DTR_ALIGN_COORD_PLOT             = '{}.dtr.aligns.png'.format(DTR_ALIGN_PREFIX)
DTR_ALIGN_TSV                    = '{}.dtr.aligns.tsv'.format(DTR_ALIGN_PREFIX)
DTR_ALIGN_CYC_PERM_TSV           = str(POLISH_DIR / 'polished.cyclic_permut.stats.tsv')

######################################
# Combine Medaka polished references #
######################################
ALL_POL_PREFIX                   = str(POLISH_DIR / 'polished')
ALL_POL_UNTRIMMED                = '{}.untrimmed.fasta'.format(ALL_POL_PREFIX)
ALL_POL_CDS_SUMMARY              = '{}.cds.summary.tsv'.format(ALL_POL_PREFIX)
ALL_POL_STRANDS                  = '{}.pol_strands.tsv'.format(ALL_POL_PREFIX)
ALL_POL_STRAND_ANNOTS            = '{}.pol_strands.reads.tsv'.format(ALL_POL_PREFIX)
ALL_POL                          = '{}.seqs.fasta'.format(ALL_POL_PREFIX)
ALL_POL_UNIQ                     = '{}.seqs.unique.fasta'.format(ALL_POL_PREFIX)
ALL_POL_DTR_STATS                = '{}.dtr.stats.tsv'.format(ALL_POL_PREFIX)
ALL_POL_DTR_GC_STATS             = '{}.dtr.gc.tsv'.format(ALL_POL_PREFIX)
ALL_POL_NUCMER_PREFIX            = '{}.nuc'.format(ALL_POL_PREFIX)
ALL_POL_NUCMER_DELTA             = '{}.nuc.delta'.format(ALL_POL_PREFIX)
ALL_POL_NUCMER_COORDS            = '{}.nuc.coords'.format(ALL_POL_PREFIX)
ALL_POL_STATS                    = '{}.stats.tsv'.format(ALL_POL_PREFIX)
ALL_POL_STATS_UNIQ               = '{}.stats.unique.tsv'.format(ALL_POL_PREFIX)
ALL_POL_CDS_PLOT_UNIQ_ALL        = '{}.unique.cds.all.png'.format(ALL_POL_PREFIX)
ALL_POL_CDS_PLOT_UNIQ_DTR_NPOL10 = '{}.unique.cds.dtr_npol10.png'.format(ALL_POL_PREFIX)

#######################################
# ANALYSIS OF LINEAR CONCATEMER READS #
#######################################
CONCATEMER_DIR                                  = OUTPUT_ROOT / SAMPLE / STYPE / VERSION / 'concatemers'
CONCATEMER_ALIGN_TMP_DIR                        = CONCATEMER_DIR / 'aln_tmp'
CONCATEMER_READ_INFO                            = str(CONCATEMER_DIR / 'concats.tsv')
CONCATEMER_READ_FASTA                           = str(CONCATEMER_DIR / 'concats.fasta')
CONCATEMER_READ_COPY_REPEATS_CONTOURS           = str(CONCATEMER_DIR / 'concats.contours.png')


wildcard_constraints:
    rank = '\d+',


##### load rules #####

include: 'rules/summary.smk'
include: 'rules/dtr_reads.smk'
include: 'rules/virsorter.smk'
include: 'rules/kaiju.smk'
include: 'rules/kmer_bins.smk'
include: 'rules/align_clusters.smk'
include: 'rules/polish.smk'
include: 'rules/annotate.smk'
include: 'rules/qc_genomes.smk'
include: 'rules/dedup.smk'
include: 'rules/dtr_align.smk'
include: 'rules/linear_concats.smk'

#############################################
# the "all" rule
#############################################
# 
# By default, all steps are run
#
output_files = [

    # all_kmer_count_and_bin: bin reads
    SUMMARY_PLOT,
    KMER_FREQS_UMAP_QSCORE,
    KMER_FREQS_UMAP_GC,
    KMER_FREQS_UMAP_READLENGTH,
    KMER_FREQS_UMAP_BINS_PLOT,

    # all_populate_kmer_bins
    lambda w: expand_template_from_bins(w, BIN_READLIST),
    lambda w: expand_template_from_bins(w, BIN_FASTA),

    # all_alignment_clusters
    KMER_BIN_STATS,
    lambda w: expand_template_from_bins(w, ALN_CLUST_OUTPUT_HEATMAP),

    # all_polish_and_annotate
    ALN_CLUST_READS_COMBO,
    lambda w: expand_template_from_bin_clusters(w, BIN_CLUSTER_REF_READ_FASTA),
    lambda w: expand_template_from_bin_clusters(w, BIN_CLUSTER_POL_READS_FASTA),
    lambda w: expand_template_from_bin_clusters(w, DTR_ALIGN_COORD_PLOT),
    lambda w: expand_template_from_bin_clusters(w, \
                                        BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT),
    lambda w: expand_template_from_bin_clusters(w, \
                                        BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS),
    lambda w: expand_template_from_bin_clusters(w, \
                                        BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS),
    lambda w: expand_template_from_bin_clusters(w, \
                                       BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS),

    # all_combine_dedup_summarize
    ALL_POL_CDS_PLOT_UNIQ_ALL,
    ALL_POL_CDS_PLOT_UNIQ_DTR_NPOL10,
    ALL_POL,
    ALL_POL_UNIQ,
    ALL_POL_STATS,

    # all_linear_concatemer_reads
    CONCATEMER_READ_COPY_REPEATS_CONTOURS,
    CONCATEMER_READ_FASTA,
]

# all_kaiju              (optional for taxonomic annotation of UMAP plots)
if config['KAIJU'].get('run', True):
    if os.path.exists(KAIJU_DB_DIR):
        output_files.extend([
            KAIJU_RESULTS_KRONA_HTML,
            expand(str(KMER_FREQS_UMAP_TAX), database=DATABASE_NAME, rank=TAX_RANK),
        ])
    else:
        logger.warning(f"No kaiju DB found in {KAIJU_DB_DIR}.\nSkipping Kaiju")


rule all:
    input:
        output_files

## The orignal partial workflows

rule all_kmer_count_and_bin:
    input:
        SUMMARY_PLOT,
        KMER_FREQS_UMAP_QSCORE,
        KMER_FREQS_UMAP_GC,
        KMER_FREQS_UMAP_READLENGTH,
        KMER_FREQS_UMAP_BINS_PLOT,

if os.path.exists(KAIJU_DB_DIR):
    rule all_kaiju:
        input:
            KAIJU_RESULTS_KRONA_HTML,
            expand(str(KMER_FREQS_UMAP_TAX), database=DATABASE_NAME, rank=TAX_RANK),
else:
    logger.warning(f"No kaiju DB found in {KAIJU_DB_DIR}.\nSkipping Kaiju")
    rule all_kaiju:
        run:
            print(f"Cannot run kaiju! DB not found at {KAIJU_DB_DIR}")

rule all_populate_kmer_bins:
    input:
        bin_reads=lambda w: expand_template_from_bins(w, BIN_READLIST),
        bin_fasta=lambda w: expand_template_from_bins(w, BIN_FASTA),

rule all_alignment_clusters:
    input:
        stats=KMER_BIN_STATS,
        heatmaps=lambda w: expand_template_from_bins(w, ALN_CLUST_OUTPUT_HEATMAP),
        align=lambda w: expand_template_from_bins(w, ALN_CLUST_OUTPUT_INFO),

rule all_polish_and_annotate:
    input:
        # all_polish_and_annotate
        ALN_CLUST_READS_COMBO,
        lambda w: expand_template_from_bin_clusters(w, BIN_CLUSTER_REF_READ_FASTA),
        lambda w: expand_template_from_bin_clusters(w, BIN_CLUSTER_POL_READS_FASTA),
        lambda w: expand_template_from_bin_clusters(w, DTR_ALIGN_COORD_PLOT),
        lambda w: expand_template_from_bin_clusters(w, \
                                            BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT),
        lambda w: expand_template_from_bin_clusters(w, \
                                            BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS),
        lambda w: expand_template_from_bin_clusters(w, \
                                            BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS),
        lambda w: expand_template_from_bin_clusters(w, \
                                           BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS),

rule all_combine_dedup_summarize:
    input:
        ALL_POL_CDS_PLOT_UNIQ_ALL,
        ALL_POL_CDS_PLOT_UNIQ_DTR_NPOL10,
        ALL_POL,
        ALL_POL_UNIQ,
        ALL_POL_STATS

rule all_linear_concatemer_reads:
    input:
        CONCATEMER_READ_COPY_REPEATS_CONTOURS,
        CONCATEMER_READ_FASTA,

