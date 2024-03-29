% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust_deltaW.R
\name{adjust_deltaW}
\alias{adjust_deltaW}
\title{Adjust inversion breakpoints based on changes in read direction}
\usage{
adjust_deltaW(
  WW_reads,
  WC_reads,
  WW_d,
  WC_d,
  inversions,
  genotype = c("het", "hom", "low")
)
}
\arguments{
\item{WW_reads}{A Granges object of reads (WW file)}

\item{WC_reads}{A GRanges object of reads (WC file)}

\item{WW_d}{A GRanges object with reads annotated with deltaW values (WW file)}

\item{WC_d}{A GRanges object with reads annotated with deltaW values (WC file)}

\item{inversions}{A GRanges object containing inversions of the specified genotype}

\item{genotype}{What genotype of inversion are we adjusting? "het", "hom", or "low", where low means posterior probability < confidence.}
}
\value{
A GRanges object of inversions with adjusted coordinates,
}
\description{
For each inversion, finds nearby peaks in deltaW values (approximating the greatest change in read direction) and makes those the new start and end coordinates.
}
\details{
Reads require some pre-processing (see process_reads()). This function widens each inversion so it contains 3x as many reads as before. This is where we will look for peaks in deltaW values to use as
breakpoints. The two highest peaks (for each inversion) become the new coordinates. If genotype is "low" or "hom", we use the peaks from the WW composite file. For "het", we take the intersection
of the intervals defined by the peaks in the two composite files, although if the result is less than 200 bp we give up and just take the peaks from the WW composite file. We keep adjusted inversions
only if they overlap the original inversion.
}
