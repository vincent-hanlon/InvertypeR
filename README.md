# InvertypeR

An R package for inversion genotyping (and discovery) using Strand-seq data. This is supplemented by an "Inversion visualization" section, at bottom, which uses PERL scripts to make ideograms linked to the UCSC Genome Browser.

For more information not provided below, or to cite this tool, please see the InvertypeR paper at https://doi.org/10.1186/s12864-021-07892-9. If you should run into trouble, please feel free to post an Issue or contact me at vhanlon [at] bccrc.ca

## Pre-processing

The main inputs to InvertypeR are:

1. Strand-seq BAM files for the individual of interest, placed together in a directory. The BAM files should have good Strand-seq quality.
2. A VCF file containing heterozygous SNVs for the individual of interest (not needed for haploids). If you don't have this, it's always best to use high-coverage NGS data of some kind to call SNVs before using InvertypeR. However, it's also possible to call SNVs directly from the Strand-seq libraries despite their shallow coverage. 
3. A BED file containing genomic intervals that you suspect may be inverted. This could be a handful of well-characterized inversions you want to genotype in a new individual, or it could be thousands of dubious inversions from the literature.
4. Appropriate priors. These aren't as hard to choose as you might think. 
5. A BED file containing regions that you expect to have poor-quality Strand-seq data. For humans, you can generally just use the one on this repo. For other species, you may want to make one.

### Strand-seq BAM files







Assuming your Strand-seq FASTQ libraries for an individual are now aligned, indexed BAM files, these scripts create two Strand-seq composite files. Poor-quality libraries must first be removed. To create the Watson-Watson (WW or WWCC) composite file, run `bash master_WWCC_composite.sh` in the directory containing the single-cell BAM files. Same goes for the Watson-Crick (WC or WCCW) composite file: run `bash master_WCCW_composite.sh`. You must first edit the header of each master script to set user-specific variables (e.g. # threads, directory containing scripts). Providing a VCF file of high-confidence heterozygous SNPs improves composite file creation, since it can be difficult to call good SNPs from sparse Strand-seq data alone. See "instructions.txt" for more details. 

In my experience, 30 or more libraries with at least 20 million unique reads in total is sufficient to make good composite files and call inversions. However, it may be possible to make good composite files with fewer reads, especially with the help of a good VCF file of heterozygous SNPs. On the other hand, having more reads will allow InvertypeR to genotype more & smaller inversions.

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

Results for the InvertypeR paper (https://doi.org/10.1186/s12864-021-07892-9) are from commit ca17a576fcbfeeb81ecd30c2c6c41ef4f1bc68cf. The newest commit produces essentially identical results, up to one or two minor shifts to inversion breakpoints per sample. 

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
