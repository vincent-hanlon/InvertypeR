# InvertypeR
Source code for all three components of the InvertypeR process:
1. Composite file creation
2. InvertypeR genotyping 
3. Inversion visualization

Composite file creation
-----------------------
Dependencies:
  - *tool (version we use)*
  - samtools (1.10)
  - freebayes (1.3.2)
  - bcftools (1.10.2)
  - R (3.5.1)
  - R package GenomicRanges (1.34.0)
  - R package [breakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) (1.5.1)
  - R package [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR) (0.99)
  - R package BSgenome.Hsapiens.UCSC.hg38 (1.4.1)

These scripts create two Strand-seq composite files, given a set of single-cell Strand-seq libraries for an individual (BAM format, indexed). Poor-quality libraries must first be removed. To create the Watson-Watson (WW or WWCC) composite file, run `bash master_WWCC_composite.sh` in the directory containing the single-cell BAM files. Same goes for the Watson-Crick (WC or WCCW) composite file: run `bash master_WCCW_composite.sh`. You must first edit the header of each master script to set user-specific variables (e.g. # threads, directory containing scripts). 

InvertypeR genotyping
-----------------------
Dependencies:
  - *tool (version we use)*
  - R (3.5.1)
  - R package invertyper (0.0.0.1)

This is done by an R package called "invertyper" that genotypes inversions in Strand-seq data. It will install a few further R packages as sub-dependencies.

The package implements a Bayesian binomial model to genotype inversions in Strand-seq data, which must be pre-processed into two composite files (WW and WC, the latter phased). In a sense, InvertypeR can also be used to discover inversions that were not already known to be present in the data. This can be done when many putative inversions are genotyped, for example, if all inversions recorded in dbVar are genotyped with an appropriate prior. This package can also adjust the start and end coordinates of inversions in a variety of ways. The directory "example_output" contains the raw InvertypeR output files (one file per run) used in the paper.

To install invertyper, you should make sure you have [devtools](https://cran.r-project.org/web/packages/devtools/index.html) installed first. Once that's done, run `devtools::install_github(repo="vincent-hanlon/InvertypeR", subdir="invertyper")`. If that doesn't work, try `devtools::install_git(url="https://github.com/vincent-hanlon/InvertypeR", ref="main", subdir="invertyper")`. If the dependency BreakpointR fails to install automatically, try installing it from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) beforehand.

For instructions on how to run the main invertyper wrapper function, load the package in R and write `?invertyper`.

Inversion visualization
-----------------------
Dependencies:
  - *tool (version we use)*
  - R (3.5.1)
  - R package dplyr (0.8.5)
  - R package gridExtra (2.3)
  - R package ggplot2 (3.3.0)
  - R package data.table (1.12.8)
  - R package psych (2.0.9)
  - ImageMagick (7.0.10-0)
  - python package img2pdf (0.4.0)
  - perl package PDF::API2 (2.038)
  - perl package LWP::UserAgent (6.49)
  
(Courtesy of Victor Guryev and Carl-Adam Mattsson)
These scripts can be found [here](https://github.com/mattssca/haploplotR), along with more detailed instructions. In brief, clone the repository, install the dependencies, put an InvertypeR output file in the `in/` directory, and put the two BreakpointR browserfiles for the WW and WC composite files (i.e. `sample_name.WW.CC.bam_reads.bed.gz` and `sample_name.WC.CW.bam_reads.bed.gz` from the composite file creation procedure) in the `in/bed_reads/` directory. Then run `bash haploplot_run.sh`. A PDF ideogram linked to a UCSC Genome Browser session will be created automatically.
