#' Overrides StrandPhaseR's function for loading sequence data
#'
#' This replaces the StrandPhaseR function by sneakily preventing it from reading from the BAM file, but instead giving it processed GAlignments results
#' A list of GAlignments objects, named by BAM file name, will need to be available to this function even though it's not passed from StrandPhaseR
#'
#' @param file Bamfile with aligned reads. (Or in InvertypeR, the filename of a bamfile)
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued. (Not needed for InvertypeR)
#' @param region If only a subset of the genomic regions should be loaded.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param filterAltAlign (A parameter used by StrandPhaseR but accepted and then ignored by InvertypeR).
#' @param min.mapq Minimum mapping quality when importing from BAM files. (not actually used by InvertypeR, same as filterAltAlign)
#' @param all_alignments A list of GAlignments objects for the Strand-seq libraries. Global so that there are no problems accessing it in parallel
#' @return a GRanges object
bamregion2GRanges_for_invertyper <- function(bamfile, bamindex, region = NULL, pairedEndReads = FALSE, min.mapq = 10, filterAltAlign = TRUE, galignmentslist = galignmentslist) {
    # Extracting the BAM name from the BAM file path + name

    galignment <- galignmentslist[[bamfile]]

    reads <- galignment_to_granges(galignment, purpose = "StrandPhaseR", paired_reads = pairedEndReads, region = region)

    return(reads)
}
