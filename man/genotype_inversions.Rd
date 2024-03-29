% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotype_inversions.R
\name{genotype_inversions}
\alias{genotype_inversions}
\title{Bayesian genotyper for inversions}
\usage{
genotype_inversions(
  WW_reads,
  WC_reads,
  regions,
  background,
  base_state,
  haploid_chromosomes = NULL,
  prior = c(0.33, 0.33, 0.33),
  haploid_prior = c(0.5, 0.5)
)
}
\arguments{
\item{WW_reads}{A GRanges object (or GAlignmentPairs in the PE case) containing reads for a WW composite file. See create_composite_files() or import_bam().}

\item{WC_reads}{A GRanges object (or GAlignmentPairs in the PE case) containing reads for a WC composite file. See create_composite_files() or import_bam().}

\item{regions}{A Granges object containing genomic intervals that are thought to be inversions.}

\item{background}{The fraction of background reads for the WW composite file. See WW_background().}

\item{base_state}{The strand state of the WW composite file: either "WW" (mostly + reads) or "CC" (mostly - reads).}

\item{haploid_chromosomes}{A vector of the names of chromosomes expected to be haploid (e.g., chrX and chrY in human males). Default NULL.}

\item{prior}{Vector of three prior weights for inversion genotypes. For example, c("ref","het","hom") = c(0.9,0.05,0.05). Default c(0.33,0.33,0.33).}

\item{haploid_prior}{Vector of two prior weights for inversions on haploid chromosomes (e.g., chrX and chrY in human males). For example, c("ref", "inv") = c(0.9,0.1). Default c(0.5,0.5).}
}
\value{
A dataframe of the regions, where each region is matched with the most probable genotype and the corresponding posterior probability, as well as some read counts.
}
\description{
Given two Strand-seq composite files (WC and WW) and a list of intervals, this computes the highest posterior probability of the possible (phased) genotypes.
}
\details{
Step-by-step: It standardizes the priors so that they sum to 1, and sets prior probabilities for errors that may be present in the data. Then it chooses the binomial probabilities for
the Bayesian model. The function counts the forward and reverse reads in each inversion, and then genotypes them in the haploid case (details omitted: it's the same as the diploid case really).
The read counts and the probabilites can be combined to calculate binomial (log) likelihoods of the possible strand states of each inversion in the two composite files (e.g. WW or WC).
Given the strand-states we expect (accounting for the errors that may be present in the data), these likelihoods can be used to compute more (log) likelihoods: this time, for the
genotypes REF, HET(0|1), HET(1|0), and HOM. We convert these into regular posterior probabilties and choose the highest one, with the associated genotype.

Note that we can phase inversions only because the WC composite file already has phased reads. This means that we know a 0|1 inversion on chr1 is on the same homolog as all 0|1 inversions
on chr1 in the sample, and that all chr1 1|0 inversions are on the other homolog. However, we don't know whether a 0|1 inversion on chr1 and a 0|1 inversion chr2 came from the same parent.
0|1 inversions are distinguished from 1|0 inversions based on the strand switch in the WC composite file ( WC -> WW or WC -> CC).
}
