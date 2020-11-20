# InvertypeR
Source code for all three components of the InvertypeR process:
1. Composite file creation
2. InvertypeR genotyping (soon to be an R package)
3. Inversion visualization (with links to the UCSC Genome Browser)

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

These scripts create two Strand-seq composite files, given a set of single-cell Strand-seq libraries for an individual (BAM format, indexed). Poor-quality libraries must first be removed. To create the Watson-Watson (WW or WWCC) composite file, run "bash master_WWCC_composite.sh" in the directory containing the single-cell BAM files. Same goes for the Watson-Crick (WC or WCCW) composite file: run "bash master_WCCW_composite.sh". Both master scripts must first be edited to set user-specific variables (e.g. # threads, directory containing scripts). 

InvertypeR genotyping
-----------------------
Dependencies:
  - *tool (version we use)*
  - R (3.5.1)
  - R package GenomicAlignments (1.18.1)

(Soon to be an R package). This series of functions is run using the invertyper() wrapper function, which takes as arguments two composite BAM files (as above), a set of intervals to genotype, priors, and a few other options (see source code for now). 


Inversion Visualization
-----------------------
(more to come)
