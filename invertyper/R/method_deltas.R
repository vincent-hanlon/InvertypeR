#' Adjust inversion breakpoints based on changes in read direction
#'
#' For each inversion, finds nearby peaks in deltaW values (approximating the greatest change in read direction) and makes those the new start and end coordinates.
#'
#' Reads require some pre-processing (see process_reads()). This function widens each inversion so it contains 3x as many reads as before. This is where we will look for peaks in deltaW values to use as 
#' breakpoints. The two highest peaks (for each inversion) become the new coordinates. If genotype is "low" or "hom", we use the peaks from the WW composite file. For "het", we take the intersection
#' of the intervals defined by the peaks in the two composite files, although if the result is less than 200 bp we give up and just take the peaks from the WW composite file. We keep adjusted inversions
#' only if they overlap the original inversion.
#'
#' @param WW_reads A Granges object of reads (WW file)
#' @param WC_reads A GRanges object of reads (WC file)
#' @param WW_d A GRanges object with reads annotated with deltaW values (WW file)
#' @param WC_d A GRanges object with reads annotated with deltaW values (WC file)
#' @param inversions A GRanges object containing inversions of the specified genotype
#' @param genotype What genotype of inversion are we adjusting? "het", "hom", or "low", where low means posterior probability < confidence.
#' @return A GRanges object of inversions with adjusted coordinates,
#'
#'
#' @export
adjust_deltaW <- function(WW_reads, WC_reads, WW_d, WC_d, inversions, genotype=c("het","hom","low")){


        #Next I want to merge and expand the inversions. The wider, expanded inversion intervals will be used to look for alternative breakpoints
        inversions <- GenomicRanges::sort(GenomicRanges::reduce(inversions))
print("initial reduce of Inversions")
print(length(inversions))
        #It's simpler to get rid of strands----otherwise, precede and follow misbehave
        GenomicRanges::strand(WW_reads) <- "*"

        #the width of inversions and hom in reads
        num_inversions <- pmax(GenomicRanges::countOverlaps(inversions, WW_reads),10)

        #Follow and precede don't work very well, especially not with empty GRanges, so I need to force zeros in the empty case
        #And also add a few annoying steps to deal with NAs appropriately.
        r <- data.frame(WW_reads)

        if(length(inversions) == 0){

                f_inversions <- 0
                p_inversions <- 0

        } else {

                f_inversions <- GenomicRanges::follow(inversions,WW_reads)
                p_inversions <- GenomicRanges::precede(inversions,WW_reads)

		if(length(inversions[is.na(f_inversions)]) > 0){
			f_inversions[is.na(f_inversions)] <- GenomicRanges::precede(inversions[is.na(f_inversions)], WW_reads)
		}
		
                if(length(inversions[is.na(p_inversions)]) > 0){
                        p_inversions[is.na(p_inversions)] <- GenomicRanges::follow(inversions[is.na(p_inversions)], WW_reads)
                }
	

        }

        #Calculating the start and end coordinates for the widened interval, and making sure they don't spill over onto the next chromosome
        #The new intervals are usually 3x as wide as the originals
        start_inversions <- ifelse(as.character(GenomicRanges::seqnames(WW_reads[pmax(1,f_inversions - num_inversions)]))==as.character(GenomicRanges::seqnames(inversions)), GenomicRanges::start(WW_reads[pmax(1,f_inversions - num_inversions)]),1)
        end_inversions <- ifelse(as.character(GenomicRanges::seqnames(WW_reads[pmin(length(WW_reads), p_inversions + num_inversions)]))==as.character(GenomicRanges::seqnames(inversions)),
                GenomicRanges::end(WW_reads[pmin(length(WW_reads), p_inversions + num_inversions)]), GenomeInfoDb::seqlengths(WW_reads)[as.character(GenomicRanges::seqnames(inversions))])

        #Assembling the new wide intervals
print("wide_inversions")
print(length(wide_inversions))
        wide_inversions <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(inversions), ranges=IRanges::IRanges(start = start_inversions, end = end_inversions))

        #Finding reads that overlap the wide intervals
	#And then some annoying checks and formatting in case the wide intervals overlap 0 reads, because the pair2frgm step removes a few
        WC_select <- IRanges::mergeByOverlaps(wide_inversions, WC_d)
	WC_missing <- wide_inversions[!IRanges::overlapsAny(wide_inversions,WC_d)]
	WC_replace <- cbind(IRanges::mergeByOverlaps(WC_missing, WC_missing, type="equal"), as.data.frame(matrix(20, ncol=6, nrow=length(WC_missing))))
        colnames(WC_replace) <- colnames(WC_select)
        WC_select <- rbind(WC_select, WC_replace)

print("1st wc_select")
print(length(WC_select))

        WW_select <- IRanges::mergeByOverlaps(wide_inversions, WW_d)
	WW_missing <- inversions[!IRanges::overlapsAny(wide_inversions,WW_d)]
	WW_replace <- cbind(IRanges::mergeByOverlaps(WW_missing, WW_missing, type="equal"), as.data.frame(matrix(20, ncol=6, nrow=length(WW_missing))))
        colnames(WW_replace) <- colnames(WW_select)
        WW_select <- rbind(WW_select, WW_replace)

print("1st ww_select")
print(length(WW_select))


        #Finding the two highest deltaW peaks per interval
        WW_coords <- do.call(rbind, as.list(IRanges::by(WW_select, WW_select[,1], function(x) boundaries(x))))
        WC_coords <- do.call(rbind, as.list(IRanges::by(WC_select, WC_select[,1], function(x) boundaries(x))))

print("WC_coords and WW_coords")
print(length(WC_coords))
print(length(WW_coords))

        if(genotype=="hom" | genotype=="low") {

                #The peaks give the new intervals
                new <-GenomicRanges::sort(GenomicRanges::GRanges(seqnames=WW_coords[,1], ranges=IRanges::IRanges(start=WW_coords[,2], width=WW_coords[,3])))

        }

        if(genotype=="het") {

                WW_ranges <- GenomicRanges::sort(IRanges::IRanges(start=WW_coords[,2], width=WW_coords[,3]))
                WC_ranges <- GenomicRanges::sort(IRanges::IRanges(start=WC_coords[,2], width=WC_coords[,3]))


                intersection <- GenomicRanges::pintersect(WW_ranges, WC_ranges, resolve.empty="start.x")

                #If the intersection is short (<200 bp) then we just take the WW interval by default, since in general the WW intervals are larger.
                intersected_ranges <- IRanges::IRanges(start=ifelse(GenomicRanges::width(intersection) < 200, GenomicRanges::start(WW_ranges), GenomicRanges::start(intersection)), width=ifelse(GenomicRanges::width(intersection) < 200,
                        GenomicRanges::width(WW_ranges),GenomicRanges::width(intersection)))

                new <-GenomicRanges::sort(GenomicRanges::GRanges(seqnames=WC_coords[,1], ranges=intersected_ranges))

        }

	
        #Removing new intervals if they don't overlap the original interval obtained by merging all overlapping inversions of the same genotype
        #A better way to do with would be to choose the two peaks so that the left one is left of GenomicRanges::end(inversions) and the right one is right of GenomicRanges::start(inversions)
        new <- new[GenomicRanges::pintersect(inversions,new)$hit]
        return(new)

}
