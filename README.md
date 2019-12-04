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
That's all. Ideally, `conda` shoud take care of all the remaining dependencies when each specific Snakemake `rule` is called.

## Usage

### Pipeline configuration
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
    - `KAIJU:db`: Path to Kaiju database to use for taxonomic annotation. Folder should contain the *.fmi, names.dmp, and nodes.dmp files. These be downloaded from the [Kaiju web server](http://kaiju.binf.ku.dk/server). nr+euk is recommended.
    - `KAIJU:switch`: Parameters to use for Kaiju annotation. These parameters should suffice in most cases.
    - `SKIP_BINS:<sample>:<stype>:<version>`: List of k-mer bins to skip for analysis. HDBSCAN assigns all unbinned reads to bin_id = -1, which should always be skipped. Occasionally downstream errors in the pipeline indicate that additional bad k-mer bins should also be skipped. 

2. The second section of `config.yml` ('Editing optional') contains software parameters that can be adjusted, but should perform sufficiently well in most cases using the default values.

### Pipeline execution
The full pipeline, starting from raw reads and ending with nanopore-polished phage genomes, can be executed in multiple steps. The first step provides summary plots for the input reads, identifies reads containing DTR sequences, creates k-mer count vectors, embeds these vectors into 2-d via [UMAP](https://github.com/lmcinnes/umap), and calls bins in the embedding via [HDBSCAN](https://github.com/lmcinnes/HDBSCAN). 
```
snakemake --use-conda -p -j <nproc> -r all_kmer_count_and_bin
```
Where `<nproc>` is the maximum number of cores for all tasks in the pipeline. The output files for this step are placed in a path according to the variables set in the `config.yml` file.
* <output_root>/<sample>/<stype>/<version>/__reads_summary__
    * `reads.summary.stats.png`: Read length and qscore distributions for all input reads
* <output_root>/<sample>/<stype>/<version>/__dtr_reads__
    * `output.dtr.fasta`: FASTA file of all DTR-containing reads
    * `output.dtr.hist.png`: Read length distribution of all DTR-containing reads
    * `output.dtr.stats.tsv`: Statistics for each DTR-containing read
* <output_root>/<sample>/<stype>/<version>/__kmer_binning__
    * `bin_membership.tsv`: bin assignments for each DTR-containing read
    * `seq_comp.*.png`: Variety of scatter plots of UMAP embedding of k-mer count vectors, annoted by features including GC-content, read length, and bin assigned by [HDBSCAN](https://github.com/lmcinnes/HDBSCAN)
    * `seq_comp.*.umap.tsv`: x- and y-coordinates for each read in the 2-d embedding
    * `seq_comp.*.umap.bins.tsv`: Mean x- and y-coordinates for each bin assigned by HDBSCAN (for finding bins in 2-d embedding)
    * `seq_comp.k5.tsv`: 5-mer count vectors for each DTR-containing read

The next step annotates reads using Kaiju and is not strictly required for producing polished genomes. However, it can be informative for verifying the integrity of the k-mer bins and for other downstream analyses.
```
snakemake --use-conda -p -j <nproc> -r all_kaiju
```
Some of the output files for this step include:
* <output_root>/<sample>/<stype>/<version>/__kaiju__
    * `results.html`: Krona dynamic plot of annotated taxonomic composition of DTR-containing reads
    * `results.taxa.tsv`: Per-read annotation results
* <output_root>/<sample>/<stype>/<version>/__kmer_binning__
    * `seq_comp.*.umap.nr_euk.[0-6].png`: Scatter plots of UMAP embedding of k-mer count vectors, annotated by various levels of taxonomic annotation

The next step simply populates subdirectories with the binned reads as assigned by HDBSCAN.
```
snakemake --use-conda -p -j <nproc> -r all_populate_kmer_bins
```
These subdirectories are located in the `kmer_binning` directory:
* <output_root>/<sample>/<stype>/<version>/kmer_binning/__bins/<bin_id>__
* Each bin subdirectory contains a list of binned read names (`read_list.txt`) and an associated FASTA file of reads (`<bin_id>.reads.fa`)

Next, each k-mer bin is refined by all-vs-all aligning all reads within a bin. The resulting alignment scores are clustered hierarchically and refined alignment clusters are called from the clustering.
```
snakemake --use-conda -p -j <nproc> -r all_alignment_clusters
```
The bin refinement results for each k-mer bin are also placed in `kmer_binning` directory:
* <output_root>/<sample>/<stype>/<version>/kmer_binning/__bins_refine/align_clusters/<bin_id>__
* Each bin refinement procedure generates an alignment (`<bin_id>.ava.paf`), clustering heatmap (`<bin_id>.clust.heatmap.png`), and information on alignment cluster assignments (`<bin_id>.clust.info.tsv`)

Next, a single read is selected from each valid alignment cluster and is polished by the remaining reads in the alignment cluster. Polishing is first done using multiple rounds of [Racon](https://github.com/isovic/racon) (3x by default), then is finished using a single round of [Medaka](https://github.com/nanoporetech/medaka) polishing.
```
snakemake --use-conda -p -j <nproc> -r all_polishing
```
This step produces polished output in the following directories:

XXX
```
snakemake --use-conda -p -j <nproc> -r all_finish_up
```
XXX

## Licence and Copyright

(c) 2019 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

## FAQs and tips
- TBD

## References and Supporting Information

A version of this pipeline is described in a preprint:

Beaulaurier J., Luo E., Eppley J., Den Uyl P., Dai X., Turner D.J., Pendelton M., Juul S., Harrington E., DeLong E.F. Assembly-free single-molecule nanopore sequencing recovers complete virus genomes from natural microbial communities. bioRxiv, doi.org/10.1101/619684 (2019).
