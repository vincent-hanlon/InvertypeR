#' Converting reads from GAlignments to GRanges formats
#'
#' The resulting GRanges objects are suitable for processing with other Strand-seq R packages like StrandPhaseR or BreakpointR
#'
#' @param galignments A GenomicAlignments object
#' @param purpose Either 'StrandPhaseR' or 'BreakpointR'. Different purposes result in different formats.
#' @param paired_reads Boolean. Are the reads paired? Default TRUE.
#' @param region A GRanges object, to which the GAlignments object can be subset.
#' @param pair2frgm Boolean. Should paired reads be merged into a single fragment, for BreakpointR? Default FALSE.
#'
#' @return A GRanges object
#'
#' @export
galignment_to_granges <- function(galignments = NULL, purpose = NULL, paired_reads = TRUE, region = NULL, pair2frgm = FALSE) {

    stopifnot("purpose must be either 'StrandPhaseR' or 'BreakpointR'." = (purpose == "StrandPhaseR" || purpose == "BreakpointR"))
    stopifnot("galignments must be of class 'GAlignmentPairs' or 'GAlignments'" = (is(galignments, "GAlignmentPairs") || is(galignments, "GAlignments")))
    stopifnot("Set paired_reads=TRUE when galignments is a GAlignmentPairs object, and set paired_reads=FALSE when galignments is a GAlignments object" = is(galignments, "GAlignmentPairs") == paired_reads)
    stopifnot("pair2frgm can only be TRUE for paired-end reads used for BreakpointR" = (pair2frgm && paired_reads && purpose == "BreakpointR") || !pair2frgm)

    if (purpose == "StrandPhaseR" && paired_reads && !pair2frgm) {
        granges.first <- GenomicRanges::granges(GenomicAlignments::first(galignments), use.mcols=TRUE)
        granges.last <- GenomicRanges::granges(GenomicAlignments::last(galignments), use.mcols=TRUE)
        granges.first$cigar <- GenomicAlignments::cigar(GenomicAlignments::first(galignments))
        granges.last$cigar <- GenomicAlignments::cigar(GenomicAlignments::last(galignments))
        GenomicRanges::strand(granges.last) <- GenomicRanges::strand(granges.first)
        granges <- GenomicRanges::sort(c(granges.first, granges.last), ignore.strand = TRUE)
        granges$XA <- NA
    } else if (purpose == "StrandPhaseR" && !paired_reads && !pair2frgm) {
        granges <- GenomicRanges::granges(galignments, use.mcols=TRUE)
        granges$cigar <- GenomicAlignments::cigar(galignments)
    } else if (purpose == "BreakpointR" && paired_reads && pair2frgm) {
        galignments.prop.pairs <- galignments[GenomicAlignments::isProperPair(galignments)]
        granges.first <- GenomicRanges::granges(GenomicAlignments::first(galignments.prop.pairs), use.mcols=TRUE)
        granges.last <- GenomicRanges::granges(GenomicAlignments::last(galignments.prop.pairs), use.mcols=TRUE)

        granges.first.plus <- granges.first[GenomicRanges::strand(granges.first) == "+"]
        granges.first.minus <- granges.first[GenomicRanges::strand(granges.first) == "-"]
        granges.last.plus <- granges.last[GenomicRanges::strand(granges.last) == "+"]
        granges.last.minus <- granges.last[GenomicRanges::strand(granges.last) == "-"]
        frag.plus.mapq <- granges.first.plus$mapq + granges.last.minus$mapq
        frag.minus.mapq <- granges.first.minus$mapq + granges.last.plus$mapq

        granges.frag.plus <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(granges.first.plus), ranges = IRanges::IRanges(start = GenomicRanges::start(granges.first.plus), end = GenomicRanges::end(granges.last.minus)), strand = GenomicRanges::strand(granges.first.plus), mapq = frag.plus.mapq)
        GenomeInfoDb::seqlengths(granges.frag.plus) <- GenomeInfoDb::seqlengths(granges.first)
        granges.frag.minus <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(granges.first.minus), ranges = IRanges::IRanges(start = GenomicRanges::start(granges.last.plus), end = GenomicRanges::end(granges.first.minus)), strand = GenomicRanges::strand(granges.first.minus), mapq = frag.minus.mapq)
        GenomeInfoDb::seqlengths(granges.frag.minus) <- GenomeInfoDb::seqlengths(granges.first)

        # Maybe this sort step could be omitted? Then it probably wouldn't take so long
        granges <- GenomicRanges::sort(c(granges.frag.plus, granges.frag.minus), ignore.strand = TRUE)
    } else if (purpose == "BreakpointR" && paired_reads && !pair2frgm) {
        granges <- as(GenomicAlignments::first(galignments), "GRanges")
        GenomicRanges::mcols(granges)[!names(GenomicRanges::mcols(granges)) %in% "mapq"] <- NULL
    } else if (purpose == "BreakpointR" && !paired_reads && !pair2frgm) {
        granges <- as(galignments, "GRanges")
        GenomicRanges::mcols(granges)[!names(GenomicRanges::mcols(granges)) %in% "mapq"] <- NULL
    }

    if (length(region) > 0 && purpose == "StrandPhaseR") {
        granges <- granges[as.vector(GenomicRanges::seqnames(granges)) %in% as.vector(GenomeInfoDb::seqlevels(region))]
        GenomeInfoDb::seqlevels(granges) <- GenomeInfoDb::seqlevels(region)
        granges <- GenomeInfoDb::keepSeqlevels(granges, GenomeInfoDb::seqlevels(region), pruning.mode = "coarse")
    }

    if (length(region) > 0 && purpose == "BreakpointR") {
        granges <- granges[S4Vectors::queryHits(GenomicRanges::findOverlaps(granges, region))]
    }

    return(granges)
}
