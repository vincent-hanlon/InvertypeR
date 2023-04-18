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

If you have all these ready, skip ahead to 'Installing and running InvertypeR', below.

### Strand-seq BAM files

Strand-seq is a single-cell library preparation for DNA sequencing that is quite good at finding inversions. Here are some relevant papers:

  - The [original](doi.org/10.1038/nmeth.2206) Strand-seq protocol
  - A recent [inversion discovery](doi.org/10.1038/s41587-020-0719-5) paper

If you have Strand-seq FASTQ files but are unsure how to align them etc., see the file instructions.txt for a step-by-step guide. 

If you have single-cell Strand-seq BAM files ready, then the next step is to select good-quality libraries. That is best done with [ASHLEYS-QC](https://github.com/friendsofstrandseq/ashleys-qc) for human libraries (for non-humans, it's probably best to have the libraries QCed manually by an expert). Install ASHLEYS-QC, and run it with default parameters as follows:

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
callvariants.sh list=samples.list ref=/PATH/TO/REFERENCE.fasta out=strand_seq_snps.vcf \
    ploidy=2 calldel=f callins=f sample="YOUR_SAMPLE_NAME"
```

InvertypeR doesn't need you to remove homozygous SNVs or indels from any VCF you provide first. It uses the package [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR) for phasing, which will figure all that out.

### Regions to genotype

InvertypeR is primarily an inversion GENOTYPING tool---users provide the coordinates of putative inversions. Although it is possible to use InvertypeR without providing any regions to genotype (in which case [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) will be used to identify putative inversions), this will tend to miss small inversions. 

If there are specific inversions you want to genotype, the choice is easy. However, if you want to find as many inversions as possible, then providing a list of coordinates from the literature or the coordinates of the spacers between pairs of inverted repeats (which commonly cause inversions) is a good idea. The file sup_table_24_inversions.bed in the example_output folder contains both of those things for GRCh38 in humans.

Regions to genotype should be provided to InvertypeR as a BED file (at least for the main function, `invertyper_pipeline()`; lower-level functions typically expect BED files to be loaded into GRanges objects with the function `import_bed()`).

### Priors

InvertypeR is a very simple Bayesian model that outputs posterior probabilties for inversion genotypes at fixed coordinates. This means we need prior probabilities. In practice, we need a vector that looks like `c(0.98, 0.01, 0.01)`, for example, which has prior probabilities that an average genomic interval in the list of regions to genotype is homozygous reference (no inversion, usually very likely), heterozygous, or homozygous alternate (both copies inverted). 

If the inversions you wish to genotype are well-characterized and you have information about allele frequencies in the population, then this is easy. This might be the case if you are only genotyping the inversions from a paper like [this one](doi.org/10.1038/s41587-020-0719-5). Use the average fraction of individuals in the study that are reference homozygotes, heterozygotes, and alternate homozygotes, respectively, for the prior vector.

If you are using the inversion list from the file sup_table_24_inversions.bed (for GRCh38 in humans), then I recommend setting `prior=c(0.9866, 0.0067, 0.0067)`
and `haploid_prior=c(0.9866, 0.0134)`. This is based on an estimate of the number of unique intervals in the inversion list. If you have your own list of putative inversions that are not well-characterized (ESPECIALLY if they're overlapping and dubious), then the InvertypeR function `choose_priors()` might help.

If you have made a list of putative inversions by calling strand switches with [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) or equivalent, then I typically use `prior=c(0.9, 0.05, 0.05)` and `haploid_prior=c(0.9, 0.1)`, but this is just based on previous experience and may differ for your particular case.

The `haploid_prior` argument is there for chrX and chrY in human males (or more generally, the sex chromosomes of the heterogametic sex), or for all chromosomes if the Strand-seq libraries are from haploid cells. In the fully haploid cells case, the arguments `chromosomes` and `haploid_chromosomes` should both be a list of all chromosomes of interest for inversion genotyping.

### Hard masked regions / blacklist

Some regions of the genome tend to have poor-quality Strand-seq data. The most problematic are things like reference assembly collapses, where only 1 copy of a region is present in the reference but several copies of the region are present in the physical genome. Depending on where the additional copies are located, the reads in the 1 reference copy map in unexpected directions. Such regions should be provided to InvertypeR as a BED file.

For GRCh38 in humans, the file blacklist.GRCh38.humans.bed is an appropriate choice. For other species, it would be best to make a hard mask file based on read depth (see the sup mat of the [InvertypeR paper](doi.org/10.1186/s12864-021-07892-9)) or based on the orientation of reads in Strand-seq libraries in many individuals.

## Installing and running InvertypeR

InvertypeR can be installed from GitHub with [devtools](https://cran.r-project.org/web/packages/devtools/index.html):
```
devtools::install_github(repo="vincent-hanlon/InvertypeR")
```

If that doesn't work, try `devtools::install_github(url="https://github.com/vincent-hanlon/InvertypeR")`. If the dependency BreakpointR fails to install automatically, try installing it from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) beforehand. Likewise, StrandPhaseR may or may need to be installed separately using `devtools::install_github()`. 

Once everything is installed and the various inputs are assembled, running InvertypeR is straightforward. The function `invertyper_pipeline()` will create composite files for an individual, attempt to discover putative inversions using [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) and the genotype them if desired, and of course genotype a user-provide list of putative inversions as well. It can also save the composite files as .RData files and write any inversions it finds to UCSC Genome Browser files along with the relevant Strand-seq reads. 

Many of these steps can also be run separately: for example, `create_composite_files()` creates composite files and `invertyper()` does the actual genotyping once composite files and putative inversions (genomic regions to genotype) are identified. 

One important argument is `adjust_method`, which controls how InvertypeR tries to adjust the provided inversion coordinates before or after genotyping. This is useful because published inversion coordinates are often inaccurate, so both inverted and non-inverted Strand-seq reads will be caught during the genotyping step. If you want to genotype the exact coordinates you provided with no alterations, choose `adjust_method="raw"`. But InvertypeR can also try to find better coordinates for putative inversions by looking for switches in the orientation of Strand-seq reads nearby the coordinates provided. This can be done only for inversions that had confident non-reference genotypes when the original coordinates were used (`adjust_method="deltas"`), or it can even be done for less confident inversions to see whether there was a problems with the coordinates provided (`adjust_method="all"`).

For more detail about the various functions in the InvertypeR package, for example for the main wrapper script, write `?invertyper_pipeline` or `?create_composite_files`, etc.

Results for the [InvertypeR paper](https://doi.org/10.1186/s12864-021-07892-9) are from commit ca17a576fcbfeeb81ecd30c2c6c41ef4f1bc68cf, when composite files were created using BASH scripts rather than the R package. The newest commit produces essentially identical results and is faster.

### Example of how to run InvertypeR

Open R and move to a directory that contains paired-end BAM files and a VCF of snps. We're assuming here that the individual is a humnan male, the BAM, VCF, and BED files all have coordinates for the reference genome GRCh38, and we want to call inversions on the sex chromosomes as well as chr1, chr2, and chr8 for some reason. 

```
library(invertyper)
i <- invertyper_pipeline(regions_to_genotype="sup_table_24_inversions.including_half_intervals.bed", 
    hard_mask='blacklist.highdepth.centromeres.bed',  prior=c(0.9866,0.0067, 0.0067), 
    haploid_prior=c(0.9866,0.0134), chromosomes=c('chr1','chr2','chr8','chrX','chrY'), 
    haploid_chromosomes=c('chrX','chrY'), vcf='snps.vcf.gz', numCPU=12, adjust_method='all', 
    save_composite_files=T, write_browser_files=T, discover_breakpointr_inversions=T,)
```

This will genotype the list of inversion provided, try to find additional inversions with the help of [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html), and save composite files and UCSC Genome Browser files.

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

These scripts can be found [here](https://github.com/mattssca/haploplotR), along with more detailed instructions. In brief, clone the repository, install the dependencies, put an InvertypeR output file in the `in/` directory, and put the browserfiles for the WW and WC composite files (from `write_browser_files()`, below) in the `in/bed_reads/` directory. Then run `bash haploplot_run.sh`. A PDF ideogram linked to a UCSC Genome Browser session will be created automatically.

The browserfiles can now be created as follows, if `create_composite_files()` or `invertyper_pipeline()` were run with `save_composite_files=T`:

```
load("./WW_composite_file.RData")
load("./WC_composite_file.RData")
write_UCSC_browser_files(WW_reads=WW_composite_file, WC_reads=WC_composite_file, paired_reads=TRUE)
```
