% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strandPhaseR_for_invertyper.R
\name{strandPhaseR_for_invertyper}
\alias{strandPhaseR_for_invertyper}
\title{A modified implementation of the core StrandPhaseR function, "strandPhaseR"}
\usage{
strandPhaseR_for_invertyper(
  numCPU = 4,
  positions = NULL,
  WCregions = NULL,
  chromosomes = NULL,
  paired_reads = TRUE,
  num.iterations = 3,
  galignmentslist = galignmentslist
)
}
\arguments{
\item{numCPU}{An integer, the number of threads to use}

\item{positions}{The path to a VCF file containing heterozygous SNPs}

\item{WCregions}{A pre-processed GRanges object that contains intervals in which specified libraries have strand-state WC (Watson-Crick)}

\item{chromosomes}{A character vector of chromosome names to be used}

\item{paired_reads}{Boolean. Are the reads paired-end?}

\item{num.iterations}{An integer. See StrandPhaseR for details. The default here is usually fine.}

\item{galignmentslist}{A list of GenomicAlignments objects containing the input Strand-seq libraries}
}
\value{
A GRanges object such that the input WCregions are annotated with 'wc' or 'cw' according to their phase
}
\description{
This function is effectively a stripped-down version of David Porubsky's strandPhaseR, but with some additional post-processing of my own to return a nice GRanges of phased WC regions.
It sneakily uses the galignmentlist (BAM files loaded into R) rather then reading them from the files all over again. This is accomplished by overwriting some StrandPhaseR functions.
}
