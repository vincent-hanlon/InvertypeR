#' Counts reads in genomic regions
#'
#' For each composite file, this function tallies up the forward and reverse reads in each putative inversion (region).
#'
#' One of the slowest steps of the genotyping process when composite files are large.
#'
#' @param WW_bam A GRanges object for the WW composite file
#' @param WC_bam A GRanges object for the WC composite file
#' @param regions A GRanges object for putative inversions
#' @return A dataframe with read counts by file
count_regions <- function(WW_reads, WC_reads, regions) {

    # Extracting counts of W and C reads by region
    WW_w <- GenomicRanges::countOverlaps(subject = WW_reads[GenomicRanges::strand(WW_reads) == "-"], query = regions)
    WW_c <- GenomicRanges::countOverlaps(subject = WW_reads[GenomicRanges::strand(WW_reads) == "+"], query = regions)
    WC_w <- GenomicRanges::countOverlaps(subject = WC_reads[GenomicRanges::strand(WC_reads) == "-"], query = regions)
    WC_c <- GenomicRanges::countOverlaps(subject = WC_reads[GenomicRanges::strand(WC_reads) == "+"], query = regions)

    # Assemmbling a data frame with region info
    counts <- data.frame(GenomicRanges::seqnames(regions), GenomicRanges::start(regions), GenomicRanges::end(regions), WW_w, WW_c, WC_w, WC_c)
    colnames(counts) <- c("chr", "start", "end", "WW_w", "WW_c", "WC_w", "WC_c")

    return(counts)
}
