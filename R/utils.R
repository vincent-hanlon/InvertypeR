
#' A utility function to import BED files
#'
#' Reads the first 3 columns of a BED file into a GRanges object.
#'
#' @param path Path to a BED file
#' @return A GRanges object
#'
#'
#' @export
import <- function(path) {

bed <- read.table(path)
gr <- GenomicRanges::sort(GenomicRanges::GRanges(seqnames=bed[,1], ranges=IRanges::IRanges(start=bed[,2], end=bed[,3])))
return(gr)

}




#' @param list_of_breakpoint_objects
#' @param states
#' @param region_size
#'
#' @return  [...]
#'
#' @export
find_regions_with_strand_state <- function(list_of_breakpoint_objects=list(), states=c('wc'), region_size=100000){

stopifnot("Allowable states are 'ww', 'cc', and 'wc'." = (!any(!states%in%c('wc','ww','cc'))))

process_counts <- function(index, states, region_size, list_of_breakpoint_objects){

	y <- list_of_breakpoint_objects[[index]]$counts
	y <- y[width(y) >= region_size & y$states %in% states]
	y$filename <- names(list_of_breakpoint_objects)[index]
	return(y)
}


regions <- lapply(c(1:length(list_of_breakpoint_objects)), process_counts, states=states, region_size=region_size, list_of_breakpoint_objects=list_of_breakpoint_objects)
regions <- sort(do.call('c', regions))

return(regions)

}





#' @param galignmentslist
#' @param regions
#' @param paired_reads
#' @param states
#' @param filename
#' @param flip_reads
#'
#' @return  [...]
extract_reads_by_region_and_state <- function(galignmentslist=NULL, regions=NULL, paired_reads=TRUE, states='wc', filename=NULL, flip_reads=FALSE) {

stopifnot("Allowable states are 'ww', 'cc', 'wc', and 'cw'." = (!any(!states%in%c('wc','ww','cc','cw'))))

reads <- galignmentslist[filename][[1]]
regions <- sort(regions[as.vector(GenomicRanges::mcols(regions)$filename)==filename & as.vector(GenomicRanges::mcols(regions)$states)==states])
reads <- reads[S4Vectors::queryHits(GenomicRanges::findOverlaps(reads, regions))]

if(flip_reads && paired_reads){
reads <- GenomicAlignments::GAlignmentPairs(GenomicAlignments::second(reads),GenomicAlignments::first(reads))
} else if(flip_reads && !paired_reads) {
reads <- BiocGenerics::invertStrand(reads)
} 

reads <- reads[!is.na(GenomicRanges::seqnames(reads))]

return(reads)

}





#' Removes the 'seq' and 'qual' mcols from a GAlignments or GAlignmentsPairs object
#'
#' I'm concerned that these fields are burdensome to store, so I want to get rid of them in the standard implementation
#'
#' @param galignments
#' @param paired_reads
#'
#' @return  [...]
#'
#' @export
drop_seq_qual <- function(galignments, paired_reads=TRUE){

stopifnot("Set paired_reads=TRUE when galignments is a GAlignmentPairs object, and set paired_reads=FALSE when galignments is a GAlignments object" = is(galignments,'GAlignmentPairs') == paired_reads)

if(!paired_reads){

	GenomicRanges::mcols(galignments)$seq <- NULL
	GenomicRanges::mcols(galignments)$qual <- NULL

} else {

	first_galignments <- GenomicAlignments::first(galignments)
	second_galignments <- GenomicAlignments::second(galignments)
	GenomicRanges::mcols(first_galignments)$seq <- NULL
	GenomicRanges::mcols(first_galignments)$qual <- NULL
        GenomicRanges::mcols(second_galignments)$seq <- NULL
        GenomicRanges::mcols(second_galignments)$qual <- NULL

	galignments <- GenomicAlignments::GAlignmentPairs(first_galignments, second_galignments)

}

return(galignments)

}





#' Finds the intersection of a set of overlapping intervals
#'
#' Used in an experimental inversion adjustment routine.
#'
#' For a GRanges object containing multiple clusters of overlapping intervals, this function acts on each cluster. For each one, it finds the largest interval that overlaps every interval in the 
#' cluster. When this is not possible, it returns a placeholder interval of length zero.
#'
#' @param interval A GRanges object
#' @return A GRanges object
#'
#'
#' @export
minimal <- function(interval){

#Group by overlaps
grouped <- GenomicRanges::findOverlaps(interval,GenomicRanges::reduce(interval),select="all")
grouped <- data.frame(S4Vectors::queryHits(grouped), S4Vectors::subjectHits(grouped))

#Separate sets of overlapping intervals and find the latest start and the earliest end
separated <- split(grouped,grouped[,2])
ranges <- cbind(sapply(c(1:length(separated)), function(x) as.character(GenomicRanges::seqnames(interval[separated[[x]][1,1]]))), sapply(c(1:length(separated)), function(x) max(GenomicRanges::start(interval[separated[[x]][,1]]))), sapply(c(1:length(separated)), function(x) min(GenomicRanges::end(interval[separated[[x]][,1]]))))

#If there is no minimal interval, replace it with an interval of length 0
ranges[as.numeric(ranges[,2]) >= as.numeric(ranges[,3]),3] <- 1
ranges[as.numeric(ranges[,2]) >= as.numeric(ranges[,3]),2] <- 1

#Creating a new set of intervals
new <- GenomicRanges::GRanges(seqnames=ranges[,1], ranges=IRanges::IRanges(start=as.numeric(ranges[,2]),end=as.numeric(ranges[,3])))

return(new)

}







#' Finds peaks in deltaW values
#'
#' Intermediate function that processes the input and output of the twopeaks function.
#'
#' Extracts a vector of deltaW values from the dataframe of GRanges objects, feeds it to twopeaks, and then extracts the genomic coordinates of the peaks that are returned.
#'
#' @param df A dataframe produced by GenomicRanges::mergeByOverlaps. Column 1 contains Granges objects (inversions), column 2 contains Granges objects that overlap the 
#'   objects in column 1 (reads) and the remaining columns are \code{mcols}, of which one must be \code{$deltaW}.
#' @return A dataframe containing genomic coordinates of the two peaks in deltaW values
#'
boundaries <- function(df){


#If there are no inversions, x$deltaW will have length zero; need to return an empty dataframe
if(length(df$deltaW)==0){

	return(data.frame(seqnames=character(), start=character(), width=character()))

}else if(length(df$deltaW)==1){

	return(data.frame(seqnames=as.character(GenomicRanges::seqnames(df[,1])),start=GenomicRanges::start(df[,1]),width=GenomicRanges::width(df[,1]) ))

}

peaks <- twopeaks(df$deltaW)

#Choosing input for IRanges based on twopeaks output
seqnames <- as.character(GenomicRanges::seqnames(df[,2][peaks[1]]))
start <- GenomicRanges::end(df[,2][peaks[1]])
end <- GenomicRanges::start(df[,2][peaks[2]])
width <- pmax(end - start + 1, 1)

return(data.frame(seqnames, start, width))

}






