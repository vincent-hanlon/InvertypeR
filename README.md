# InvertypeR
v1.0.0.1

InvertypeR is an R package for genotyping (and discovering) inversions using Strand-seq data. This is supplemented by an "Inversion visualization" section, at bottom, which uses R/PERL scripts to make ideograms linked to the UCSC Genome Browser.

For further details beyond what is provided below, or to cite InvertypeR, please see the [InvertypeR paper](https://doi.org/10.1186/s12864-021-07892-9). If you should run into trouble, please feel free to post an Issue or contact me at vincent [AT] alumni.ubc.ca 
## Quick start guide

Install InvertypeR from GitHub using the R package devtools:
```
devtools::install_github(repo="vincent-hanlon/InvertypeR")
```

Then collect good-quality Strand-seq libraries for a single diploid individual in an input directory (see below for haploid species), along with a VCF file of their SNVs. 
Provide a list of putative inversions (for humans: hanlon_2021_BMCgenomics_augmented.bed from this repo) and appropriate priors. Provide a list of hard masked regions (for humans: hard_mask.GRCh38.humans.bed from this repo).

Run the full InvertypeR pipeline using something like this command, which will probably take several hours:
```
library(invertyper)
i <- invertyper_pipeline(regions_to_genotype="hanlon_2021_BMCgenomics_augmented.bed", 
    hard_mask='hard_mask.GRCh38.humans.bed', prior=c(0.9683, 0.0158, 0.0158), 
    haploid_prior=c(0.9683, 0.0317), chromosomes=c('chr1','chr2','chr8','chrX','chrY'), 
    haploid_chromosomes=c('chrX','chrY'), vcf='snps.vcf.gz', numCPU=12, adjust_method='all', 
    save_composite_files=T, write_browser_files=T, discover_breakpointr_inversions=T,)
```

If your diploid individual has haploid sex chromosomes, list them under the haploid_chromosomes argument as well as under chromosomes (for human males: `haploid_chromosomes=c("chrX", "chrY")`)

InvertypeR outputs a tab-delimited text file with the posterior probabilties for inversion genotypes. Optionally, you can also save the composite files to an RData file (since they are slow to make), and you can write reads and inversion calls to UCSC genome browser files to examine the data yourself. InvertypeR also outputs PDF files containing [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) plots of the composite files, so that you can check they generated correctly (described below).

## Long form user guide

- [Pre-processing](#pre-processing)
  * [Strand-seq BAM files](#strand-seq-bam-files)
  * [VCF file](#vcf-file)
  * [Regions to genotype](#regions-to-genotype)
  * [Priors](#priors)
  * [Hard masked regions / blacklist](#hard-masked-regions--blacklist)
- [Running InvertypeR](#running-invertyper)
  * [Installation and overview](#installation-and-overview)
  * [Example of how to run InvertypeR](#example-of-how-to-run-invertyper)
  * [Soft masked regions: checking your composite files](#soft-masked-regions-checking-your-composite-files)
  * [Interpreting InvertypeR output](#interpreting-invertyper-output)
  * [A note on ploidy](#a-note-on-ploidy)
- [Inversion visualization](#inversion-visualization)

### Pre-processing

The main inputs to InvertypeR are:

1. Strand-seq BAM files for the individual of interest, placed together in a directory. The BAM files should have good Strand-seq quality.
2. A VCF file containing heterozygous SNVs for the individual of interest (not needed for haploids). If you don't have this, it is usually possible to call enough SNVs from the Strand-seq data.
3. A BED file containing genomic intervals that you suspect may be inverted. This could be a handful of well-characterized inversions you want to genotype in a new individual, or it could be thousands of dubious inversions from the literature.
4. Appropriate priors. These aren't as hard to choose as you might think. 
5. A BED file containing regions that you expect to have poor-quality Strand-seq data. For humans, you can generally just use the one on this repo. For other species, you may want to make one.

If you have all these ready, you can probably skip ahead to 'Installing and running InvertypeR', below.

#### Strand-seq BAM files

Strand-seq is a single-cell library preparation for DNA sequencing that is quite good at finding inversions. Here are some relevant papers:

  - The [original](https://doi.org/10.1038/nmeth.2206) Strand-seq protocol
  - A recent [inversion discovery](https://doi.org/10.1038/s41587-020-0719-5) paper

If you have Strand-seq FASTQ files but are unsure how to align them etc., see the file alignment_and_QC_instructions.txt for a step-by-step guide. 

If you have single-cell Strand-seq BAM files ready, then the next step is to select good-quality libraries. That is best done with [ASHLEYS-QC](https://github.com/friendsofstrandseq/ashleys-qc) for human libraries aligned to the GRCh38 reference genome (for non-humans, it's probably best to have the libraries QCed manually by an expert). Install ASHLEYS-QC, and run it with default parameters as follows:

```
/PATH/TO/FILE/ashleys.py -j 12 features -f ./ -w 5000000 2000000 1000000 800000 600000 400000 200000 -o ./features.tsv
/PATH/TO/FILE/ashleys.py predict -p ./features.tsv -o ./quality.txt -m /PATH/TO/OTHER/FILE/svc_default.pkl
```

Then, libraries that score lower than 0.5 are considered poor quality and should be moved to a subdirectory. 

Typically, I use at least 30 libraries with at least 20 million non-duplicate aligned reads in total for inversion genotyping. However, having more reads will allow InvertypeR to genotype more (smaller) inversions.

#### VCF file

InvertypeR uses two composite files consisting of reads merged from all the Strand-seq libraries. Phase/haplotype information is required to create one of the two composite files: luckily, Strand-seq is very good at phasing. However, Strand-seq libraries have shallow depth of coverage, so they are not ideal for calling high-confidence heterozygous SNVs. If WGS sequence data is available for the individual you wish to genotype, it is best to use that to call SNVs and provide the VCF file as input to InvertypeR. 

However, if that is not possible, then it is usually feasible to call enough SNVs from both good- and poor-quality Strand-seq libraries using call_variants.sh from [bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/), which can be installed using [conda](https://anaconda.org/agbiome/bbtools). The resulting SNVs can contain lots of false positives, although I have generally been able to make good composite files anyway.

To install:
```
conda install -c agbiome bbtools=37*
```

To run:
```
ls ./*.bam bam_poor_quality/*.bam > samples.list
callvariants.sh list=samples.list ref=/PATH/TO/REFERENCE.fasta out=strand_seq_snps.vcf \
    ploidy=2 calldel=f callins=f sample="YOUR_SAMPLE_NAME"
```

InvertypeR doesn't need you to remove homozygous SNVs or indels from any VCF you provide first. It uses the package [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR) for phasing, which will figure all that out.

#### Regions to genotype

InvertypeR is primarily an inversion GENOTYPING tool---users provide the coordinates of putative inversions. Although it is possible to use InvertypeR without providing any regions to genotype (in which case [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) will be used to identify putative inversions), this will tend to miss small inversions. 

If there are specific inversions you want to genotype, the choice is easy. However, if you want to find as many inversions as possible, then providing a list of coordinates from the literature or the coordinates of the spacers between pairs of inverted repeats (which commonly cause inversions) is a good idea. The file hanlon_2021_BMCgenomics_augmented.bed contains both of those things for GRCh38 in humans. Alternatively, if you want to genotype a smaller set of common inversions in human populations, the 292 inversions in the file porubsky_2022_cell.bed are a good place to start, which are from [this paper](https://doi.org/10.1016/j.cell.2022.04.017). See below to choose appropriate priors.

Regions to genotype should be provided to InvertypeR as a BED file (at least for the main function, `invertyper_pipeline()`; lower-level functions typically expect BED files to be loaded into GRanges objects with the function `import_bed()`).

#### Priors

InvertypeR is a very simple Bayesian model that outputs posterior probabilties for inversion genotypes at fixed coordinates. This means we need prior probabilities. In practice, we need a vector that looks like `c(0.98, 0.01, 0.01)`, for example, which has prior probabilities that an average genomic interval in the list of regions to genotype is homozygous reference (no inversion, usually very likely), heterozygous, or homozygous alternate (both copies inverted). 

If the inversions you wish to genotype are well-characterized and you have information about allele frequencies in the population, then this is easy. This might be the case if you are only genotyping the inversions from a paper like [this one](https://doi.org/10.1016/j.cell.2022.04.017). Use known allele frequencies, or the average fraction of individuals in the study that are reference homozygotes, heterozygotes, and alternate homozygotes, respectively.

The `haploid_prior` argument is there for chrX and chrY in human males (or more generally, the sex chromosomes of the heterogametic sex), or for all chromosomes if the Strand-seq libraries are from haploid cells.

<strong>Table 1. </strong><em>Some recommended priors for provided regions to genotype (for humans, GRCh38). The first two rows were done using InvertypeR's `choose_priors()` function, the third row was done using the genotypes of simple inversions from [Porubsky et al. 2022](https://doi.org/10.1016/j.cell.2022.04.017), and the last row is a rule of thumb used used for the arguments `breakpointr_prior` and `breakpointr_haploid_prior`.</em>
| `regions_to_genotype`  | `prior` | `haploid_prior` |
| ------------- | ------------- | ------------- |
| hanlon_2021_BMCgenomics_augmented.bed |  `c(0.9683, 0.0158, 0.0158)` | `c(0.9683, 0.0317)` |
| hanlon_2021_BMCgenomics_original.bed |  `c(0.9730, 0.0135, 0.0135)` | `c(0.9730, 0.0270)` |
| porubsky_2022_cell.bed  | `c(0.5166, 0.4043, 0.0791)` | `c(0.6202, 0.3798)` |
| BreakpointR as used by InvertypeR | `c(0.9, 0.05, 0.05)` | `c(0.9, 0.1)` | 

#### Hard masked regions / blacklist

Some regions of the genome tend to have poor-quality Strand-seq data. The most problematic are things like reference assembly collapses, where only 1 copy of a region is present in the reference but several copies of the region are present in the physical genome. Depending on where the additional copies are located, the reads in the 1 reference copy map in unexpected directions. Such regions should be provided to InvertypeR as a BED file (for `invertyper_pipeline()`; lower-level functions expext a GRanges object produced with `import_bed()`).

For GRCh38 in humans, the file hard_mask.GRCh38.humans.bed is an appropriate choice. For other species, it would be best to make a hard mask file based on read depth (see the sup mat of the [InvertypeR paper](https://doi.org/10.1186/s12864-021-07892-9)) or based on the orientation of reads in Strand-seq libraries in many individuals.

### Running InvertypeR

Results for the [InvertypeR paper](https://doi.org/10.1186/s12864-021-07892-9) are from commit ca17a576fcbfeeb81ecd30c2c6c41ef4f1bc68cf, when composite files were created using BASH scripts rather than the R package. The newest commit (described in this user guide) produces essentially identical results and is faster and simpler.

#### Installation and overview

InvertypeR can be installed from GitHub with [devtools](https://cran.r-project.org/web/packages/devtools/index.html):
```
devtools::install_github(repo="vincent-hanlon/InvertypeR")
```

If that doesn't work, try `devtools::install_github(url="https://github.com/vincent-hanlon/InvertypeR")`. If the dependency BreakpointR fails to install automatically, try installing it from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) beforehand. Likewise, StrandPhaseR may or may need to be installed separately using `devtools::install_github()`. 

Once everything is installed and the various inputs are assembled, running InvertypeR is straightforward. The function `invertyper_pipeline()` will create composite files for an individual, attempt to discover putative inversions using [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) and the genotype them if desired, and of course genotype a user-provide list of putative inversions as well. It can also save the composite files as .RData files and write any inversions it finds to UCSC Genome Browser files along with the relevant Strand-seq reads. 

Many of these steps can also be run separately: for example, `create_composite_files()` creates composite files and `invertyper()` does the actual genotyping once composite files and putative inversions (genomic regions to genotype) are identified. 

One important argument is `adjust_method`, which controls how InvertypeR tries to adjust the provided inversion coordinates before or after genotyping. This is useful because published inversion coordinates are often inaccurate, so both inverted and non-inverted Strand-seq reads will be caught during the genotyping step. If you want to genotype the exact coordinates you provided with no alterations, choose `adjust_method="raw"`. But InvertypeR can also try to find better coordinates for putative inversions by looking for switches in the orientation of Strand-seq reads nearby the coordinates provided. This can be done only for inversions that had confident non-reference genotypes when the original coordinates were used (`adjust_method="deltas"`), or it can even be done for less confident inversions to see whether there was a problems with the coordinates provided (`adjust_method="all"`).

For more detail about the various functions in the InvertypeR package, for example for the main wrapper script, write `?invertyper_pipeline` or `?create_composite_files`, etc.

#### Example of how to run InvertypeR

Open R and move to a directory that contains paired-end BAM files and a VCF of snps. We're assuming here that the individual is a human male, the BAM, VCF, and BED files all have coordinates for the reference genome GRCh38, and we want to call inversions on the sex chromosomes as well as chr1, chr2, and chr8 for some reason. 

```
library(invertyper)
i <- invertyper_pipeline(regions_to_genotype="hanlon_2021_BMCgenomics_augmented.bed", 
    hard_mask='hard_mask.GRCh38.humans.bed',  prior=c(0.9683, 0.0158, 0.0158), 
    haploid_prior=c(0.9683, 0.0317), chromosomes=c('chr1','chr2','chr8','chrX','chrY'), 
    haploid_chromosomes=c('chrX','chrY'), vcf='snps.vcf.gz', numCPU=12, adjust_method='all', 
    save_composite_files=T, write_browser_files=T, discover_breakpointr_inversions=T,)
```

This will genotype the list of inversion provided, try to find additional inversions with the help of [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html), and save composite files and UCSC Genome Browser files.

#### Soft masked regions: checking your composite files

Composite file creation isn't perfect. The function `create_composite_files()` generates a PDF of each composite file created, which you should inspect for anomalies. Most commonly, a large (real) inversion or reference misorient can reorient so many reads that InvertypeR has trouble identifying strand states nearby, so it just doesn't use reads from that region or chromosome. The solution is to find the coordinates of those regions (typically using [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html)), and provide them to InvertypeR using the `soft_mask` argument. Reads from those regions will appear in composite files and be used for inversion calling (so the large inversion/misorient won't be overlooked), but they won't be used to identify strand states to MAKE the composite files. This often appears as a chromosome that has way fewer reads than other chromosomes (see below). However, some variation in coverage between chromosomes is expected, especially when a low-coverage chromosome in one composite file (e.g., Watson-Watson, or WW) has high coverage in the other composite file (e.g., Watson-Crick or WC). 

For humans, the big inversion on chr8 very occasionally causes this problem, or in some males the large inversion at the start of chrY. The coordinates for these using GRCh38 are chr8:8129404-12608295 and chrY:6240762-9910296, respectively. Even if you don't expect to have composite file creation issues, it would be reasonable to put them in a BED file and pass them to InvertypeR using `soft_mask` as a preventative measure. 

<img src="https://github.com/vincent-hanlon/InvertypeR/blob/main/composite_errors.png" width=680>
<strong>Figure 1. </strong><em>An example of a non-human WW composite file with an error caused by reference assembly misorients on chromosome 16 (green bumps). </em>
<br>
<br>

Other rare issues: If you know that a large heterozygous inversion (>2 Mb) is present in your sample but it doesn't appear in the WC composite file (or it appears, but with an unclear strand switch that has high background), something may be wrong with phasing. Alternatively, very occasionally the WW composite file might have high background, that is, across a large region or chromosome, there is a consistent fraction of reads which are not the same orientation (or colour) as the majority. 

Finally, it is generally a good idea to check that your two composite files have roughly equal numbers of reads, and that combined they have a similar (though lower) number of reads than individual Strand-seq BAMs you provided as input. When counting reads, be sure to exclude duplicate and unmapped reads, and to distinguish counts of read pairs from counts of individual reads.

#### Interpreting InvertypeR output

The first 3 columns of InvertypeR's primary output file (below) are the genomic coordinates of the inversions. The next 4 columns (WW_counts_C, etc.) are counts of forward and reverse reads in the two composite files, i.e., the information InvertypeR used to make genotype calls. Column 8 is the inversion genotype, including reference aka non-inverted genotypes. Heterozygous inversions on the same chromosome are phased relative to each other. Column 9 is the posterior probability of the genotype (and typically we only consider a genotype 'confident' if the posterior is above some threshold, like 0.95), and column 10 indicates whether the inversion coordinates span a large read-poor area like a centromere, and may need to be adjusted manually based using the UCSC browser files.

<br>

<img src="https://github.com/vincent-hanlon/InvertypeR/blob/main/output_example.png">
<strong>Figure 2. </strong><em>InvertypeR's primary output file.</em>

<br>

#### A note on ploidy

InvertypeR was built for diploids, and although it won't work well on polyploids, haploids are ok. The difference here is that we don't need a Watson-Crick composite file or a VCF of SNVs to build one with. 

The arguments `chromosomes` and `haploid_chromosomes` are a little convoluted in this case. For haploids, they should be equal---they should both be a list of all chromosomes of interest for inversion genotyping. If you're running `invertyper_pipeline()` that should be all you need to know, but if you're running `create_composite_files()` you need to set `type='ww'`.

### Inversion visualization

Beyond the UCSC Genome Browser files and [BreakpointR](https://bioconductor.org/packages/release/bioc/html/breakpointR.html) plots generated by the InvertypeR R package itself, it is possible to create chromosome ideograms displaying the location, genotype, and phase of all called inversions as clickable links that display the relevant region and reads in the genome browser. These links usually expire after a week or so, but they can be quite useful for exploring inversion calls. The linked ideograms are created by running accessory R/bash/perl scripts as described below.

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

The required browserfiles can now be created as follows, if `create_composite_files()` or `invertyper_pipeline()` were run with `save_composite_files=T`:

```
load("./WW_composite_file.RData")
load("./WC_composite_file.RData")
write_UCSC_browser_files(WW_reads=WW_composite_file, WC_reads=WC_composite_file, paired_reads=TRUE)
```
