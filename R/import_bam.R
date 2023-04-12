#' Reads a BAM file into a GAlignments object, subsetting to regions/chromosomes of interest
#'
#' It can also take a GAlignments or GRanges object, in which case it will just apply the necessary region subsetting instead of trying to load it fresh somehow.
#'
#' @param bam Path to a BAM file. Alternatively, a GRanges or GAlignments object
#' @param region A GRanges object of regions from which reads should be taken. Default NULL.
#' @param paired_reads Boolean. Are the reads paired?
#' @param blacklist A GRanges object of regions from which reads should NOT be taken. Default NULL
#' @param chromosomes A character vector of chromosome names. Reads from other chromosomes will not be loaded. Default NULL (all chromosomes).
#'
#' @return A GAlignments object of reads
#'
#' @export
import_bam <- function(bam = NULL, region = NULL, paired_reads = TRUE, blacklist = NULL, chromosomes = NULL) {

    if (length(region) == 0 & is.character(bam)) {
        lengths <- Rsamtools::scanBamHeader(bam)[[1]]$targets
        region <- GenomicRanges::GRanges(names(lengths), ranges = IRanges::IRanges(1, lengths))
    } else if (length(region) > 0) {
        region <- GenomicRanges::reduce(region, min.gapwidth = 1000)
    } else if (length(region) == 0 & !is.character(bam) & (!is.null(chromosomes) | !is.null(blacklist))) {
        lengths <- GenomeInfoDb::seqlengths(bam)
        region <- GenomicRanges::GRanges(seqnames = names(lengths), ranges = IRanges::IRanges(start = c(1:length(lengths)), end = lengths))
    } else {
        warning("import_bam is returning the input reads unaltered because neither a subsetting region nor the path to a BAM file were supplied. This may also be a consequence of not providing a blacklist")
        return(bam)
    }

    if (length(blacklist) > 0) {
        region <- GenomicRanges::setdiff(region, blacklist)
    }

    if (!is.null(chromosomes)) {
        region <- region[as.vector(GenomicRanges::seqnames(region)) %in% chromosomes]
    }

    if (is.character(bam)) {
        stopifnot("Set paired_reads=FALSE when your BAM file contains single-end reads, and vice versa" = suppressMessages(Rsamtools::testPairedEndBam(bam)) == paired_reads)

        if (paired_reads) {
            # Reading in files. Warnings suppressed because GenomicAlignments::readGAlignment... generates them for no good reason
            # Note that the intervals were already extended by 1Mb (so that we can read in flanking region for other purposes) and reduced so they aren't overlapping
            p1 <- Rsamtools::ScanBamParam(
                flag = Rsamtools::scanBamFlag(isPaired = TRUE, isUnmappedQuery = FALSE, isDuplicate = FALSE, isSecondaryAlignment = FALSE, hasUnmappedMate = FALSE, isProperPair = TRUE),
                mapqFilter = 10, which = region, what = c("seq", "qual", "mapq", "cigar", "flag")
            )
            reads <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bam, param = p1))
        } else {
            p1 <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired = FALSE, isUnmappedQuery = FALSE, isDuplicate = FALSE, isSecondaryAlignment = FALSE), mapqFilter = 10, which = region, what = c("seq", "qual", "mapq", "cigar", "flag"))
            reads <- suppressWarnings(GenomicAlignments::readGAlignments(bam, param = p1))
        }
    } else if (length(region) > 0) {
        reads <- bam[S4Vectors::queryHits(GenomicAlignments::findOverlaps(bam, region))]
    } else {
        stop("The import_bam function requires either the path to a BAM file, or a GenomicAlignments object full of reads AND a GRanges object containing genomic regions for which the corresponding reads will be selected.")
    }

    return(reads)
}
