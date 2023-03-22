
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
#' @return A list of four GRanges objects: the input objects WW_reads and WC_reads and the two output objects (for WW and WC composite files) that have been annotated with deltaW values.
#'
#'
#' @export
process_reads <- function(WW_reads, WC_reads, sex, paired_reads=TRUE){

        #Processing input reads for Aaron's deltaWCalculator

        if (paired_reads) {

                WW_reads <- galignment_to_granges(WW_reads, purpose='BreakpointR', paired_reads=TRUE, pair2frgm=TRUE)
                WC_reads <- galignment_to_granges(WC_reads, purpose='BreakpointR', paired_reads=TRUE, pair2frgm=TRUE)


        } else {

                WW_reads <- GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(WW_reads),ranges=IRanges::IRanges(start=GenomicRanges::start(WW_reads), end=GenomicRanges::end(WW_reads)), strand=GenomicRanges::strand(WW_reads), GenomicRanges::mcols(WW_reads), 
seqlengths=GenomeInfoDb::seqlengths(WW_reads))
                WC_reads <- GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(WC_reads),ranges=IRanges::IRanges(start=GenomicRanges::start(WC_reads), end=GenomicRanges::end(WC_reads)), strand=GenomicRanges::strand(WC_reads), GenomicRanges::mcols(WC_reads), 
seqlengths=GenomeInfoDb::seqlengths(WC_reads))

        }

	# in rare cases when a few reads in a female misalign to chrY, deltaWCalculator will sometimes return objects of different lengths for different binsizes, which causes a warning
	# this avoids the issue
	if (sex == "female") {

		WW_reads <- WW_reads[GenomicRanges::seqnames(WW_reads)!="chrY"]
                WC_reads <- WC_reads[GenomicRanges::seqnames(WC_reads)!="chrY"]

	}


        #Calculating deltaW values
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

        #Taking an arbitrary sum of deltaWs to look near and far (including 2x the deltaW10)
        WW_d$deltaW <- WW_d$deltaW + 2*WW_d_10$deltaW + WW_d_20$deltaW + WW_d_40$deltaW + WW_d_80$deltaW
        WC_d$deltaW <- WC_d$deltaW + 2*WC_d_10$deltaW + WC_d_20$deltaW + WC_d_40$deltaW + WC_d_80$deltaW



        return(list(WW_reads, WC_reads,WW_d,WC_d))

}



