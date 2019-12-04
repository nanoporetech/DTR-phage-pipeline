from pathlib import Path
import os
from glob import glob

configfile: 'config.yml'
SCRIPT_DIR = srcdir('scripts')

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
OUTPUT_ROOT  = Path(config['output_root'])
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
SUMMARY_PLOT   = str(READS_DIR / 'reads.summary.stats.png')


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
# KAIJU CLASSIFICATION            #
###################################
KAIJU_DIR                 = OUTPUT_ROOT /  SAMPLE / STYPE / VERSION / 'kaiju'
KAIJU_RESULTS             = str(KAIJU_DIR / 'results.tsv')
KAIJU_RESULTS_TAXA        = str(KAIJU_DIR / 'results.taxa.tsv')
KAIJU_RESULTS_KRONA       = str(KAIJU_DIR / 'results.krona')
KAIJU_RESULTS_KRONA_HTML  = str(KAIJU_DIR / 'results.html')
KAIJU_TAXIDS              = str(KAIJU_DIR / 'taxids.csv')
KAIJU_TAXINFO             = str(KAIJU_DIR / 'taxinfo.csv')
KAIJU_TAX_GENOMESIZE_EST  = str(KAIJU_DIR / 'taxid_genomesize_estimates.csv')


###################################
# K-mer UAMP files and plots      #
###################################
KMER_BIN_ROOT                     = OUTPUT_ROOT /  SAMPLE / STYPE / VERSION / 'kmer_binning'
KMER_BINS_MEMBERSHIP              = str(KMER_BIN_ROOT / 'bin_membership.tsv')
KMER_BIN_STATS                    = str(KMER_BIN_ROOT / 'bin_stats.csv')
KMER_FREQS_TMP                    = str(KMER_BIN_ROOT / 'seq_comp.k{}.tmp'.format(K))
KMER_FREQS                        = str(KMER_BIN_ROOT / 'seq_comp.k{}.tsv'.format(K))
KMER_FREQS_UMAP                   = str(KMER_BIN_ROOT / 'seq_comp.k{}.{}bp.umap.tsv'.format(K,config['UMAP']['min_rl']))
KMER_FREQS_UMAP_TAX               = str(KMER_BIN_ROOT / ('seq_comp.k%s.%sbp.umap.{database}.{rank}.png' % (K,config['UMAP']['min_rl'])))
KMER_FREQS_UMAP_QSCORE            = str(KMER_BIN_ROOT / 'seq_comp.k{}.{}bp.umap.qscore.png'.format(K,config['UMAP']['min_rl']))
KMER_FREQS_UMAP_GC                = str(KMER_BIN_ROOT / 'seq_comp.k{}.{}bp.umap.gc.png'.format(K,config['UMAP']['min_rl']))
KMER_FREQS_UMAP_READLENGTH        = str(KMER_BIN_ROOT / 'seq_comp.k{}.{}bp.umap.readlength.png'.format(K,config['UMAP']['min_rl']))
KMER_FREQS_UMAP_BINS_PLOT         = str(KMER_BIN_ROOT / 'seq_comp.k{}.{}bp.umap.bins.png'.format(K,config['UMAP']['min_rl']))
KMER_FREQS_UMAP_BINS_COORDS       = str(KMER_BIN_ROOT / 'seq_comp.k{}.{}bp.umap.bins.tsv'.format(K,config['UMAP']['min_rl']))


###################################
# BIN ANALYSIS FILES              #
###################################
BINS_DIR                    = KMER_BIN_ROOT / 'bins'
BIN_READLIST                = str(BINS_DIR / '{bin_id}' / 'read_list.txt')
BIN_FASTA                   = str(BINS_DIR / '{bin_id}' / '{bin_id}.reads.fa')
BIN_GENOME_SIZE             = str(BINS_DIR / '{bin_id}' / 'genome_size.txt')
ALL_BIN_IDS                 = list(map(lambda x: x.split('/')[-1], glob(str(BINS_DIR)+'/*')))
SKIP_BINS                   = config['SKIP_BINS'][SAMPLE][STYPE][VERSION]
BIN_IDS                     = [b for b in ALL_BIN_IDS if int(b) not in SKIP_BINS]
BINNED_ANALYSIS_ROOT        = KMER_BIN_ROOT / "bins_refine"

#####################################
# Alignment clustering  FILES       #
#####################################
ALN_CLUST_DIR            = BINNED_ANALYSIS_ROOT / 'align_clusters'
ALN_CLUST_PAF            = str(ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.ava.paf')
ALN_CLUST_OUTPUT_PREFIX  = str(ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust')
ALN_CLUST_OUTPUT_HEATMAP = str(ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust.heatmap.png')
ALN_CLUST_OUTPUT_INFO    = str(ALN_CLUST_DIR / '{bin_id}' / '{bin_id}.clust.info.csv')
ALN_CLUST_READS_COMBO    = str(ALN_CLUST_DIR / 'all_bins.clust.info.csv')


###################################
# Separate cluster-specific reads #
###################################
BIN_CLUSTER_DIR                 = BINNED_ANALYSIS_ROOT / 'bin_cluster_reads'
BIN_CLUSTER_READS_INFO          = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.readinfo.csv')
BIN_CLUSTER_READS_LIST          = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.readlist.csv')
BIN_CLUSTER_READS_FASTA         = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.reads.fasta')
BIN_CLUSTER_REF_READ_LIST       = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_readlist.csv')
BIN_CLUSTER_POL_READS_LIST      = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.pol_readlist.csv')
BIN_CLUSTER_REF_READ_FASTA      = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.fa')
BIN_CLUSTER_POL_READS_FASTQ     = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.pol_reads.fq')
BIN_CLUSTER_POL_READS_FASTA     = str(BIN_CLUSTER_DIR / '{bin_clust_id}' / '{bin_clust_id}.pol_reads.fa')


####################################
# Polish cluster reads using Racon #
####################################
POLISH_DIR                       = BINNED_ANALYSIS_ROOT / 'bin_cluster_polishing'
RACON_DIR                        = POLISH_DIR / 'racon'
BIN_CLUSTER_RACON_POLISHED_FASTA = str(RACON_DIR / '{bin_clust_id}' / '{{bin_clust_id}}.ref_read.racon_{repeats}x.fasta'.format(repeats=RACON_ROUNDS))


#####################################
# Polish cluster reads using Medaka #
#####################################
MEDAKA_DIR                                    = POLISH_DIR / 'medaka'
BIN_CLUSTER_POLISHED_REF_TMP                  = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.tmp.fasta')
BIN_CLUSTER_POLISHED_REF                      = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.fasta')
BIN_CLUSTER_POLISHED_POL_VS_REF_PAF           = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.paf')
BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS       = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.strands.tsv')
BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.strands.annots.tsv')
BIN_CLUSTER_POLISHED_REF_PRODIGAL             = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.fasta')
BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT         = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.txt')
BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS       = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read.medaka.prodigal.cds.stats.txt')


#########################################
# Check for fixed or cyclic permutation #
#########################################
DTR_ALIGN_PREFIX                 = str(MEDAKA_DIR / '{bin_clust_id}' / '{bin_clust_id}.ref_read')
DTR_ALIGN_HIST_PLOT              = '{}.dtr.hist.png'.format(DTR_ALIGN_PREFIX)
DTR_ALIGN_COORD_PLOT             = '{}.dtr.aligns.png'.format(DTR_ALIGN_PREFIX)
DTR_ALIGN_TSV                    = '{}.dtr.aligns.tsv'.format(DTR_ALIGN_PREFIX)
DTR_ALIGN_CYC_PERM_TSV           = str(POLISH_DIR / 'all_ref_reads.medaka.cyclic_permut.stats.tsv')


######################################
# Combine Medaka polished references #
######################################
ALL_POLISHED_PREFIX                              = str(POLISH_DIR / 'all_ref_reads.medaka')
ALL_POLISHED_COMBINED                            = '{}.fasta'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_CDS_SUMMARY                = '{}.prodigal.cds.summary.tsv'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_CDS_PLOT_ALL               = '{}.prodigal.cds.summary.all.png'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_CDS_PLOT_DTR               = '{}.prodigal.cds.summary.dtr.png'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_CDS_PLOT_DTR_NPOL10        = '{}.prodigal.cds.summary.dtr_npol10.png'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_CDS_PLOT_DTR_NPOL20        = '{}.prodigal.cds.summary.dtr_npol20.png'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_STRANDS                    = '{}.pol_strands.tsv'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_STRAND_ANNOTS              = '{}.pol_strands.annots.tsv'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_STRAND_UMAP                = '{}.pol_strands.annots.umap.png'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED                    = '{}.trimmed.fasta'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED_DTR_STATS          = '{}.trimmed.dtr.tsv'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED_DTR_GC_ABUND_STATS = '{}.trimmed.dtr.gc.abund.tsv'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED_DTR                = '{}.trimmed.dtr.fasta'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED_DTR_UNIQ           = '{}.trimmed.dtr.unique.fasta'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_PREFIX  = '{}.trimmed.dtr.nuc'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_DELTA   = '{}.trimmed.dtr.nuc.delta'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_COMBINED_TRIMMED_DTR_NUCMER_COORDS  = '{}.trimmed.dtr.nuc.coords'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_MGA_CDS_SUMMARY                     = '{}.mga.cds.summary.tsv'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_STATS                               = '{}.stats.tsv'.format(ALL_POLISHED_PREFIX)
ALL_POLISHED_STATS_FILT                          = '{}.stats.filtered.tsv'.format(ALL_POLISHED_PREFIX)


#######################################
# ANALYSIS OF LINEAR CONCATEMER READS #
#######################################
CONCATEMER_DIR                                  = OUTPUT_ROOT / SAMPLE / STYPE / VERSION / 'concatemers'
CONCATEMER_ALIGN_TMP_DIR                        = CONCATEMER_DIR / 'aln_tmp'
CONCATEMER_READ_INFO                            = str(CONCATEMER_DIR / 'concats_{}_{}.gt{}.tsv'.format(SAMPLE,STYPE,config['LIN_CONCAT']['minlen']))
CONCATEMER_READ_FASTA                           = str(CONCATEMER_DIR / 'concats_{}_{}.gt{}.fasta'.format(SAMPLE,STYPE,config['LIN_CONCAT']['minlen']))
CONCATEMER_READ_PDF_PLOT_COPIES                 = str(CONCATEMER_DIR / 'concats_{}_{}.gt{}.kde.copies.png'.format(SAMPLE,STYPE,config['LIN_CONCAT']['minlen']))
CONCATEMER_READ_PDF_PLOT_LENGTHS                = str(CONCATEMER_DIR / 'concats_{}_{}.gt{}.kde.lengths.png'.format(SAMPLE,STYPE,config['LIN_CONCAT']['minlen']))
CONCATEMER_READ_COPY_REPEATS_CONTOURS           = str(CONCATEMER_DIR / 'concats_{}_{}.gt{}.contours.png'.format(SAMPLE,STYPE,config['LIN_CONCAT']['minlen']))


wildcard_constraints:
    rank = '\d+',


##### load rules #####

include: 'rules/summary.smk'
include: 'rules/dtr_reads.smk'
include: 'rules/kaiju.smk'
include: 'rules/kmer_bins.smk'
include: 'rules/align_clusters.smk'
include: 'rules/polish.smk'
include: 'rules/prodigal.smk'
include: 'rules/qc_genomes.smk'
include: 'rules/dedup.smk'
include: 'rules/dtr_align.smk'
include: 'rules/linear_concats.smk'

#############################################
# Required steps for assembly-free analysis #
#############################################

# 1. all_kmer_count_and_bin
# 2. all_kaiju              (optional for taxonomic annotation of UMAP plots)
# 3. all_populate_kmer_bins
# 4. all_alignment_clusters
# 5. all_polishing
# 6. all_finish_up

rule all_kmer_count_and_bin:
    input:
        SUMMARY_PLOT,
        KMER_FREQS_UMAP_QSCORE,
        KMER_FREQS_UMAP_GC,
        KMER_FREQS_UMAP_READLENGTH,
        KMER_FREQS_UMAP_BINS_PLOT,

rule all_kaiju:
    input:
        KAIJU_RESULTS_KRONA_HTML,
        expand(KMER_FREQS_UMAP_TAX, database=DATABASE_NAME, rank=TAX_RANK),

rule all_populate_kmer_bins:
    input: dynamic(BIN_READLIST),
           dynamic(BIN_FASTA),

rule all_alignment_clusters:
    input:
        KMER_BIN_STATS,
        expand(ALN_CLUST_OUTPUT_HEATMAP, bin_id=BIN_IDS),
        expand(ALN_CLUST_OUTPUT_INFO, bin_id=BIN_IDS),

rule all_polishing:
    input:
        ALN_CLUST_READS_COMBO,
        dynamic(BIN_CLUSTER_READS_INFO),
        dynamic(BIN_CLUSTER_READS_LIST),
        dynamic(BIN_CLUSTER_REF_READ_LIST),
        dynamic(BIN_CLUSTER_REF_READ_FASTA),
        dynamic(BIN_CLUSTER_POL_READS_LIST),
        dynamic(BIN_CLUSTER_POL_READS_FASTA),
        dynamic(BIN_CLUSTER_RACON_POLISHED_FASTA),
        dynamic(BIN_CLUSTER_POLISHED_REF),

rule all_finish_up:
    input:
        dynamic(DTR_ALIGN_COORD_PLOT),
        dynamic(BIN_CLUSTER_POLISHED_REF_PRODIGAL),
        dynamic(BIN_CLUSTER_POLISHED_REF_PRODIGAL_TXT),
        dynamic(BIN_CLUSTER_POLISHED_REF_PRODIGAL_STATS),
        dynamic(BIN_CLUSTER_POLISHED_POL_VS_REF_STRANDS),
        dynamic(BIN_CLUSTER_POLISHED_POL_VS_REF_STRAND_ANNOTS),
        ALL_POLISHED_COMBINED_STRANDS,
        ALL_POLISHED_COMBINED_STRAND_ANNOTS,
        ALL_POLISHED_COMBINED_STRAND_UMAP,
        ALL_POLISHED_COMBINED_CDS_PLOT_ALL,
        ALL_POLISHED_COMBINED_CDS_PLOT_DTR,
        ALL_POLISHED_COMBINED_CDS_PLOT_DTR_NPOL10,
        ALL_POLISHED_COMBINED_CDS_PLOT_DTR_NPOL20,
        ALL_POLISHED_COMBINED_CDS_SUMMARY,
        ALL_POLISHED_COMBINED_TRIMMED,
        ALL_POLISHED_COMBINED_TRIMMED_DTR_STATS,
        ALL_POLISHED_COMBINED_TRIMMED_DTR_GC_ABUND_STATS,
        ALL_POLISHED_COMBINED_TRIMMED_DTR,
        ALL_POLISHED_COMBINED_TRIMMED_DTR_UNIQ,
        ALL_POLISHED_STATS,
        ALL_POLISHED_STATS_FILT,
        DTR_ALIGN_CYC_PERM_TSV,
        CONCATEMER_READ_PDF_PLOT_COPIES,
        CONCATEMER_READ_PDF_PLOT_LENGTHS,
        CONCATEMER_READ_COPY_REPEATS_CONTOURS,
        CONCATEMER_READ_FASTA,
