# InvertypeR

InvertypeR is an R package for genotyping (and discovering) inversions using Strand-seq data. This is supplemented by an "Inversion visualization" section, at bottom, which uses R/PERL scripts to make ideograms linked to the UCSC Genome Browser.

For further details beyond what is provided below, or to cite InvertypeR, please see the [InvertypeR paper](doi.org/10.1186/s12864-021-07892-9). If you should run into trouble, please feel free to post an Issue or contact me at vincent [AT] alumni.ubc.ca 

## Pre-processing

The main inputs to InvertypeR are:

1. Strand-seq BAM files for the individual of interest, placed together in a directory. The BAM files should have good Strand-seq quality.
2. A VCF file containing heterozygous SNVs for the individual of interest (not needed for haploids). If you don't have this, it is usually possible to call enough SNVs from the Strand-seq data.
3. A BED file containing genomic intervals that you suspect may be inverted. This could be a handful of well-characterized inversions you want to genotype in a new individual, or it could be thousands of dubious inversions from the literature.
4. Appropriate priors. These aren't as hard to choose as you might think. 
5. A BED file containing regions that you expect to have poor-quality Strand-seq data. For humans, you can generally just use the one on this repo. For other species, you may want to make one.

### Strand-seq BAM files

Strand-seq is a single-cell library preparation for DNA sequencing that is quite good at finding inversions. Here are some relevant papers:

  - The [original](doi.org/10.1038/nmeth.2206) Strand-seq protocol
  - A recent [inversion discovery](doi.org/10.1038/s41587-020-0719-5) paper

If you have Strand-seq FASTQ files but are unsure how to align them etc., see the file instructions.txt for a step-by-step guide. 

If you have single-cell Strand-seq BAM files ready, then the next step is to select good-quality libraries. That is best done with [ASHLEYS-QC](https://github.com/friendsofstrandseq/ashleys-qc). Install ASHLEYS-QC, and run it with default parameters as follows:

```
/PATH/TO/FILE/ashleys.py -j 12 features -f ./ -w 5000000 2000000 1000000 800000 600000 400000 200000 -o ./features.tsv
/PATH/TO/FILE/ashleys.py predict -p ./features.tsv -o ./quality.txt -m /PATH/TO/OTHER/FILE/svc_default.pkl
```

Then, libraries that score lower than 0.5 are considered poor quality and should be moved to a subdirectory. 

Typically, I use at least 30 libraries with at least 20 million non-duplicate aligned reads in total for inversion genotyping. However, having more reads will allow InvertypeR to genotype more (smaller) inversions.

### VCF file

InvertypeR uses two composite files consisting of reads merged from all the Strand-seq libraries. Phase/haplotype information is required to create one of the two composite files: luckily, Strand-seq is very good at phasing. However, Strand-seq libraries have shallow depth of coverage, so they are not ideal for calling high-confidence heterozygous SNVs. If WGS sequence data is available for the individual you wish to genotype, it is best to use that to call SNVs and provide the VCF file as input to InvertypeR. 

However, if that is not possible, then it is usually feasible to call enough SNVs from the good- and poor-quality Strand-seq libraries using [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/), which can be installed using [conda](https://anaconda.org/bioconda/bbmap):

```
ls ./*.bam bam_poor_quality/*.bam > samples.list
callvariants.sh list=samples.list ref=/PATH/TO/FILE/GRCh38.fasta out=strand_seq_snps.vcf ploidy=2 calldel=f callins=f sample="YOUR_SAMPLE_NAME"
```

InvertypeR doesn't need you to remove homozygous SNVs or indels from any VCF you provide first. It uses the package [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR) for phasing, which will figure all that out.

### Regions to genotype

InvertypeR is primarily an inversion GENOTYPING tool---users provide the coordinates of putative inversions, 

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

## Inversion visualization

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
