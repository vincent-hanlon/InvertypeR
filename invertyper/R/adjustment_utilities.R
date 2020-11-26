#' Finds two peaks in a vector
#' 
#' Looks for the two biggest peaks (local maxima) in a vector of non-negative integers.
#'
#' This is used to determine the start and end coordinates of an inversion from a list of deltaW values, which are associated with reads. It uses run length encoding 
#' to find the two highest local maxima. It doesn't select two adjacent peaks: instead, it introduces a spacer equal to 1/12th the length of the vector. It returns the
#' coordinates of the right side of the leftmost peak and the left side of the rightmost peak. 
#'
#' @param x A vector of non-negative integers.
#' @return A vector of length two.
#' @examples
#' 
#' twopeaks(c(1,2,4,2,1,3,3,3,5,5,5,6,6,2,2,2,5,4,5,4,7))
#'
#'
#'
#' @export
twopeaks <- function(x){

t <- rle(x)

if(length(t$length) > 1){

	#retain only local maxima in t
	t$values[-(which(diff(sign(diff(t$values)))==-2)+1)] <- 0

	#Finding the highest peak
	m1 <- which.max(t$values)
	
	#Setting a rough spacer: peaks shouldn't be too close together!
	space <- as.integer(length(t$values)/12)
		
	#Set highest peak and nearby values to zero
	t$values[max(m1-space,1):min(m1+space,length(t$values))] <- 0

	#Finding the second-highest peak, and figuring out which one comes first
	m2 <- which.max(t$values)
	m <- sort(c(m1,m2))

	#Returning the inner coordinates of the peaks
	coords <- c(sum(t$lengths[1:m[1]]), sum(t$lengths[1:(m[2]-1)])+1)

} else if (length(t$length) == 1){

	coords <- (c(1, length(x))) 
}

return(coords)

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
#' @export
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
grouped <- data.frame(queryHits(grouped), subjectHits(grouped))

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
