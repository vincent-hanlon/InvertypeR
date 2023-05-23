#' Read processing for inversion adjustment
#'
#' Input reads are sent through the deltaWCalculator from BreakpointR, which associates them with deltaW values approximating the local change in read direction
#'
#' For each read, we sum multiple deltaW values calculated for binsizes 5, 10, 10, 20, 40, and 80. Yes, 10 is there twice: it means we're prioritizing changes in read direction from
#' the smaller scale.
#'
#' @param WW_reads A GRanges object from the WW composite file.
#' @param WC_reads A GRanges object from the WC composite file.
#' @param paired_reads Boolean: are the reads paired-end? Default TRUE.
#' @param haploid_chromosomes A vector of the names of chromosomes expected to be haploid (e.g., chrX and chrY in human males).
#' @return A list of four GRanges objects: the input objects WW_reads and WC_reads and the two output objects (for WW and WC composite files) that have been annotated with deltaW values.
process_reads <- function(WW_reads, WC_reads, haploid_chromosomes, paired_reads = TRUE) {

    # Processing input reads for Aaron's deltaWCalculator

    WW_reads <- galignment_to_granges(WW_reads, purpose = "BreakpointR", paired_reads = paired_reads, pair2frgm = paired_reads)
    WC_reads <- galignment_to_granges(WC_reads, purpose = "BreakpointR", paired_reads = paired_reads, pair2frgm = paired_reads)

    # in rare cases when a few reads misalign to a haploid sex chromosome (e.g., when reads from human females appear on chrY), deltaWCalculator will sometimes return objects of different lengths for different binsizes, which causes a warning
    # this avoids the issue

    absent_haploid_chromosomes <- c("chrY", "chrW")[!c("chrY", "chrW") %in% haploid_chromosomes]
    WW_reads <- WW_reads[!as.vector(GenomicRanges::seqnames(WW_reads)) %in% absent_haploid_chromosomes]
    WC_reads <- WC_reads[!as.vector(GenomicRanges::seqnames(WC_reads)) %in% absent_haploid_chromosomes]

    # Calculating deltaW values
    WW_d <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 5))
    WC_d <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 5))

    WW_d_10 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 10))
    WC_d_10 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 10))

    WW_d_20 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 20))
    WC_d_20 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 20))

    WW_d_40 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 40))
    WC_d_40 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 40))

    WW_d_80 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 80))
    WC_d_80 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 80))

    # Taking an arbitrary sum of deltaWs to look near and far (including 2x the deltaW10)
    t <- tryCatch({
        WW_d$deltaW <- WW_d$deltaW + 2 * WW_d_10$deltaW + WW_d_20$deltaW + WW_d_40$deltaW + WW_d_80$deltaW
        WC_d$deltaW <- WC_d$deltaW + 2 * WC_d_10$deltaW + WC_d_20$deltaW + WC_d_40$deltaW + WC_d_80$deltaW
    }, error = function(e) {
        message("error combining deltaW values in InvertypeR's process_reads() function")
        message(e)
    }, warning = function(w) {
        message("\nWarning while combining deltaW values in invertyper::process_reads(). Did you provide the wrong haploid sex chromosomes? (e.g., setting haploid_chromosomes='chrY' for a human female)")
        message(paste0(w, "\n"))
    })

    return(list(WW_reads, WC_reads, WW_d, WC_d))
}
