![ONT_logo](/ONT_logo.png)

-----------------------------

DTR phage pipeline
==================
The DTR phage pipeline is a Snakemake pipeline for obtaining polished, full-length  dsDNA bacteriophage genomes with direct terminal repeats (DTRs) from enviromental samples. 

Getting started
===============
## Dependencies
The pipeline relies heavily  on `conda` to supply [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and the software packages required for various steps in the analysis. If `conda` is not currently installed on your system (highly encouraged!), please first follow the installation instructions for installing [miniconda](https://docs.conda.io/en/latest/miniconda.html).

[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) is a workflow management system that uses Python-based language that is ideally-suited for providing reproducible and scalable bioinformatic pipelines. If you are unfamiliar with Snakemake, this [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) is a good place to start. 

## Installation
First clone the repositorty to the desired location and enter the cloned directory:
```
git clone https://github.com/nanoporetech/DTR-phage-pipeline.git
cd DTR-phage-pipeline
```
Next, create a `conda` environment where you will run `Snakemake`:
```
conda env create -f environment.yml
```
That's all. Ideally, `conda` shoud take care of all the remaining dependencies when each specific Snakemake step (described below) is executed.

### If conda is slow

For speed tips related to conda, see *mamba* section below.

### Testing

To verify that everything is set up correctly, simply run the workflow without
modifying config.yaml:

```
snakemake --use-conda -p -r -j 10
```

## Pipeline configuration

The pipeline determines the input files, output path, and the various software parameters from the file `config.yml`, which contains key:value pairs that control how the pipeline runs. 
1. Attention must be paid to the first section of `config.yml` ('Editing mandatory'), as these might change from run to run.
    - `sample`: name of the sample (e.g. 25m)
    - `stype`: sample type (e.g. 1D or filtered_02um, etc)
    - `version`: run version (for changing run parameters)
    - `input_fasta`: FASTA file containing nanopore reads to be analyzed
    - `input_summary`: Sequencing summary file associated with the FASTA file
    - `output_root`: Base directory where pipeline results will be written
    - `max_threads`: Maximum number of threads to allocate for each rule
    - `MEDAKA:model`: Medaka model to use for polishing the draft genome (e.g. `r941_min_high_g303`). See Medaka [documentation](https://github.com/nanoporetech/medaka#models) for details on selecting the proper model for you basecalled reads.
    - `MEDAKA:threads`: Number of threads to use for each Medaka genome polishing. Should be << `max_threads` so that multiple genomes can be polished in parallel.
    - `KAIJU:db`: Path to Kaiju database to use for taxonomic annotation. Folder should contain the `*.fmi`, `names.dmp`, and `nodes.dmp` files. These be downloaded from the [Kaiju web server](http://kaiju.binf.ku.dk/server). nr+euk is recommended.
    - `KAIJU:switch`: Parameters to use for Kaiju annotation. These parameters should suffice in most cases.
    - `SKIP_BINS:<sample>:<stype>:<version>`: List of k-mer bins to skip for analysis. HDBSCAN assigns all unbinned reads to bin_id = -1, which should always be skipped. Occasionally downstream errors in the pipeline indicate that additional bad k-mer bins should also be skipped. 

2. The second section of `config.yml` ('Editing optional') contains software parameters that can be adjusted, but should perform sufficiently well in most cases using the default values.

## Pipeline execution
The full pipeline, starting from raw reads and ending with nanopore-polished phage genomes, can be executed in one step. 

```
snakemake --use-conda -p -j <nproc> -r 
```
Where `<nproc>` is the maximum number of cores for all tasks in the pipeline. The output files for this workflow are placed in a path according to the variables set in the `config.yml` file. In this description, `<run_output_dir>` will refer to the path consisting of `<output_root>/<sample>/<stype>/<version>` as defined in the `config.yml` file.

### all_kmer_count_and_bin
The first step provides summary plots for the input reads, identifies reads containing DTR sequences, creates k-mer count vectors, embeds these vectors into 2-d via [UMAP](https://github.com/lmcinnes/umap), and calls bins in the embedding via [HDBSCAN](https://github.com/lmcinnes/HDBSCAN). 

The outputs from this step will fall into three directories:

 * `<run_output_dir>/reads_summary`
    * `reads.summary.stats.png`: Read length and qscore distributions for all input reads
 * `<run_output_dir>/dtr_reads`
    * `output.dtr.fasta`: FASTA file of all DTR-containing reads
    * `output.dtr.hist.png`: Read length distribution of all DTR-containing reads
    * `output.dtr.stats.tsv`: Statistics for each DTR-containing read
 * `<run_output_dir>/kmer_binning`
    * `bin_membership.tsv`: bin assignments for each DTR-containing read
    * `kmer_comp.tsv`: 5-mer count vectors for each DTR-containing read
    * `kmer_comp.umap.tsv`: x- and y-coordinates for each read in the 2-d embedding
    * `kmer_comp.umap.*.png`: Variety of scatter plots of [UMAP](https://github.com/lmcinnes/umap) embedding of k-mer count vectors, annoted by features including GC-content, read length, and bin assigned by [HDBSCAN](https://github.com/lmcinnes/HDBSCAN)
    * `kmer_comp.umap.bins.tsv`: Mean x- and y-coordinates for each bin assigned by HDBSCAN (for finding bins in 2-d embedding)
    
### all_kaiju
The next (optional) step annotates reads using Kaiju and is not strictly required for producing polished genomes. However, it can be informative for verifying the integrity of the k-mer bins and for other downstream analyses.

Some of the output files for this step include:
* `<run_output_dir>/kaiju`
    * `results.html`: Krona dynamic plot of annotated taxonomic composition of DTR-containing reads
    * `results.taxa.tsv`: Per-read annotation results
* `<run_output_dir>/kmer_binning`
    * `kmer_comp.umap.nr_euk.[0-6].png`: Scatter plots of UMAP embedding of k-mer count vectors, annotated by various levels of taxonomic annotation

### all_populate_kmer_bins
The next step simply populates subdirectories with the binned reads as assigned by HDBSCAN.

These subdirectories are located in the `kmer_binning` directory:
* `<run_output_dir>/kmer_binning/bins/<bin_id>`
* Each bin subdirectory contains a list of binned read names (`read_list.txt`) and an associated FASTA file of reads (`<bin_id>.reads.fa`)

### all_alignment_clusters
Next, each k-mer bin is refined by all-vs-all aligning all reads within a bin. The resulting alignment scores are clustered hierarchically and refined alignment clusters are called from the clustering.

The bin refinement results for each k-mer bin are also placed in `kmer_binning` directory:
* `<run_output_dir>/kmer_binning/refine_bins/align_clusters/<bin_id>`
* Each bin refinement procedure generates an alignment (`<bin_id>.ava.paf`), clustering heatmap (`<bin_id>.clust.heatmap.png`), and information on alignment cluster assignments (`<bin_id>.clust.info.tsv`)

### all_polish_and_annotate
Next, a single read is selected from each valid alignment cluster and is polished by the remaining reads in the alignment cluster. Polishing is first done using multiple rounds of [Racon](https://github.com/isovic/racon) (3x by default), then is finished using a single round of [Medaka](https://github.com/nanoporetech/medaka) polishing.

This step produces polished output in the following directories:
* `<run_output_dir>/kmer_binning/refine_bins/align_cluster_reads/<clust_id>` simply contains the reads corresponding to each alignment cluster
* `<run_output_dir>/kmer_binning/refine_bins/align_cluster_polishing/racon/<clust_id>` contains the Racon polishing output for each alignment cluster
* `<run_output_dir>/kmer_binning/refine_bins/align_cluster_polishing/medaka/<clust_id>` contains the Medaka polishing output for each alignment cluster
    * The critical output file in each of these `medaka/<clust_id>` folders is the `<clust_id>.ref_read.medaka.fasta` file containing the __Medaka-polished genome produced from this alignment cluster__. Subsequent Snakemake rules analyze, aggregate, and deduplicate these polished genomes.
    * `<clust_id>.ref_read.medaka.prodigal.cds.*` files describe the coding sequence annotations from [Prodigal](https://github.com/hyattpd/Prodigal) for this Medaka-polished genome.
    * `<clust_id>.ref_read.strands.*` files describe the strand abundance for reads in each alignment cluster. Clusters containing >80% reads from a single strand should be treated with suspicion.
    * `<clust_id>.ref_read.dtr.aligns.*` files describe the results of aligning the DTR sequence from each corresponding k-mer bin to the polished genome from each alignment cluster. If the DTR sequences all align to the same reference positions, the DTR is fixed. However, if they align all over the reference genome, this suggests that a headful DNA packaging mechanism was used. 

### all_combine_dedup_summarize
Next, we finish up the genome discovery portion of the pipeline by running a series of aggregations and evaluations of the final polished genome sequences. 

All output from this step is written to a single directory:
* `<run_output_dir>/kmer_binning/refine_bins/align_cluster_polishing`
    * `polished.seqs.fasta`: The combined set of genomes from each alignment cluster
    * `polished.seqs.unique.fasta`: Same as above but only containing unique sequences after the deduplication step
    * `polished.stats.tsv`: Various statistics for each polished genome, including length, GC content, DTR details, cluster strand abundance, CDS annotation statistics, circular permutation status, and many others
    * `polished.stats.unique.tsv`: Same as above but only containing unique sequences after the deduplication step
    * `polished.unique.cds.summary.all.png`: Summary plots of summary statistics for the coding sequences (CDS) annotated by Prodigal for each unique polished genome
    * `polished.unique.cds.summary.dtr_npol10.png`: Same as above but only including polished genomes with a confirmed DTR and at least 10 reads used for polishing

### all_linear_concatemer_reads
Finally, we run one final step to query the sequencing dataset for linear concatemer reads that could represent interesting mobile elements in the environmental sample.

All output from this step is written to a single directory:
* `<run_output_dir>/concatemers`
    * `concats.fasta`: All identified concatemeric reads found in the input reads
    * `concats.tsv`: Statistic for each concatemeric read found, including readname, length, repeat size, and repeat copy count
    * `concats.contours.png`: Scatter plot with density contours showing the relationship between the observed repeat length and copy count in all identified concatemeric reads.

## Licence and Copyright

(c) 2019 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

## FAQs and tips

### *mamba*
The conda package manager is a powerful tool for replicating bioinformatics,
but it can be slow for large envirnonments. If you find installing any of the 
environments is taking too long, try using
[mamba](https://github.com/TheSnakePit/mamba).

Mamba is a fast drop-in replacement for conda. It is still in development, but
so far has proven to e useful. Simply install in your base envirnoment:

```
conda install mamba
```

and then replace `conda` with `mamba` in commands that install packages. EG:
```
mamba env create -f environment.yml
```

To get snakemake to use mamba, pass `--conda-frontend mamba` to your command:
```
snakemake --use-conda --conda-frontend mamba -p -r -j <nprocs>
```

### more to come
 - TBD

## References and Supporting Information

This pipeline is described in Genome Research:

Beaulaurier J., Luo E., Eppley J., Den Uyl P., Dai X., Turner D.J., Pendelton M., Juul S., Harrington E., DeLong E.F. Assembly-free single-molecule sequencing recovers complete virus genomes from natural microbial communities. Genome Research (2020). [doi:10.1101/gr.251686.119](https://genome.cshlp.org/content/early/2020/02/19/gr.251686.119.abstract)
