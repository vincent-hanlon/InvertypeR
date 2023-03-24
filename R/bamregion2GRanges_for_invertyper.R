#' Overrides StrandPhaseR's function for loading sequence data
#'
#' This replaces the StrandPhaseR function by sneakily preventing it from reading from the BAM file, but instead giving it processed GAlignments results
#' A list of GAlignments objects, named by BAM file name, will need to be available to this function even though it's not passed from StrandPhaseR
#'
#' @param bamfile
#' @param bamindex
#' @param region
#' @param pairedEndReads
#' @param min.mapq
#' @param filterAltAlign
#' @param all_alignments
#' @return a GRanges object
#' @export
bamregion2GRanges_for_invertyper <- function(bamfile, bamindex, region=NULL, pairedEndReads=FALSE, min.mapq=10, filterAltAlign=TRUE, all_alignments=galignmentslist_global_for_invertyper) {

	# Extracting the BAM name from the BAM file path + name
	message('using Vincents altered bamregion2GRanges...')

	bamfile <- intToUtf8(rev(utf8ToInt(bamfile)))
	bamfile <- gsub("/.*","", bamfile)

	bamfile <- intToUtf8(rev(utf8ToInt(bamfile)))

	galignment <- all_alignments[[bamfile]]

	reads <- galignment_to_granges(galignment, purpose="StrandPhaseR", paired_reads=pairedEndReads, region=region)
	
	name <- paste0(as.character(bamfile),".txt")
	file.create(name)

	return(reads)

}
