#' Wrapper for methods to adjust inversion breakpoints
#'
#' Provides a few options as to how to adjust inversion start and end coordinates.  
#'
#' See invertyper() for more details.
#'
#' @param inversions A dataframe of inversions output by genotype_inversions.
#' @param reads A list of two Granges objects: reads from a WW composite file and a WC composite file.
#' @param confidence Posterior probability threshold above which you consider genotype calls to be reliable. Used to decide whether to keep adjusted inversions, as well as to identify low-confidence calls 
#'   for adjust_method "low". Default 0.95.
#' @param base A list output by WWCC_background().
#' @param sex Sex of sample to figure out sex chromosomes. Default "female".
#' @param prior Vector of three prior weights for inversion genotypes. For example, c("ref","het","hom") = c(0.9,0.05,0.05).
#' @param prior_male Vector of two prior weights for male sex chromosomes. For example, c("ref", "inv") = c(0.9,0.1). Default c(0.5,0.5). 
#' @param adjust_method. See invertyper() for more details. Default "deltas". 
#' @param paired_reads Boolean: are the reads paired-end? Default TRUE.
#' @return A dataframe with adjusted inversions.
#'
#'
#' @export
adjust_intervals <- function(inversions, reads, confidence=0.95, base,  sex="female", prior, prior_male=c(0.5,0.5), adjust_method="deltas",paired_reads=TRUE){ 

		#separating called inversions by genotype and probability and creating GRanges objects; I don't bother with reference homozygotes
		het <- inversions[inversions$probability>=confidence & (inversions$genotype == "1|0" | inversions$genotype == "0|1"),]
		hom <- inversions[inversions$probability>=confidence & (inversions$genotype == "1" | inversions$genotype == "1|1") ,]
		het <- GenomicRanges::GRanges(seqnames=het[,1], ranges=IRanges::IRanges(start=het[,2], end=het[,3]), mcols= het[,-c(1:3)])
                hom <- GenomicRanges::GRanges(seqnames=hom[,1], ranges=IRanges::IRanges(start=hom[,2], end=hom[,3]), mcols=hom[,-c(1:3)])

		if (adjust_method == "deltas" | adjust_method == "all" | adjust_method=="low") {
        	
			new_reads <- process_reads(reads[[1]],reads[[2]], paired_reads=paired_reads, sex=sex)
		}

		if (adjust_method == "deltas" | adjust_method == "all") {


		        new_het <- adjust_deltaW(new_reads[[1]], new_reads[[2]], new_reads[[3]], new_reads[[4]], inversions=het, genotype="het")
			new_hom <- adjust_deltaW(new_reads[[1]], new_reads[[2]], new_reads[[3]], new_reads[[4]], inversions=hom, genotype="hom")

			delta_inversions <- genotype_new_regions(list(new_het,new_hom), het, hom, reads, confidence, base, sex, prior, prior_male)		
			
			#Combining inversions whether or not the delta method worked. Either way it's worth trying to merge those that overlap.
			het <- delta_inversions[[1]]
			hom <- delta_inversions[[2]]				
	
		}

		#Need an empty dataframe for the last rbind, in case adjust_method isn't "all"
                new_low_ref <- data.frame(character(),integer(),integer(),integer(),integer(),integer(),integer(),character(),double())
                colnames(new_low_ref) <- c("chr","start","end",paste0(base[[2]],"_counts_W"), paste0(base[[2]],"_counts_C"), "WC_counts_W", "WC_counts_C","genotype","probability")


                if (adjust_method == "all" | adjust_method == "low" ) {

                        low <- inversions[inversions$probability<confidence,]
                        low <- GenomicRanges::GRanges(seqnames=low[,1], ranges=IRanges::IRanges(start=low[,2], end=low[,3]), mcols= low[,-c(1:3)])
                        new_low <- adjust_deltaW(new_reads[[1]], new_reads[[2]], new_reads[[3]], new_reads[[4]], inversions=low, genotype="low")
                        new_low <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=new_low, background=base[[1]], base_state=base[[2]],  sex=sex, prior=prior, prior_male=prior_male)
                        new_low_het <- new_low[new_low$probability>=confidence & (new_low$genotype=="0|1" | new_low$genotype=="1|0"),]
                        new_low_hom <- new_low[new_low$probability>=confidence & (new_low$genotype=="1|1" | new_low$genotype=="1"),]
                        new_low_ref <- new_low[new_low$probability>=confidence & (new_low$genotype=="0|0" | new_low$genotype=="0"),]
		
	                new_low_het <- GenomicRanges::GRanges(seqnames=new_low_het[,1], ranges=IRanges::IRanges(start=new_low_het[,2], end=new_low_het[,3]), mcols= new_low_het[,-c(1:3)])
			new_low_hom <- GenomicRanges::GRanges(seqnames=new_low_hom[,1], ranges=IRanges::IRanges(start=new_low_hom[,2], end=new_low_hom[,3]), mcols= new_low_hom[,-c(1:3)])
			new_low_ref <- GenomicRanges::GRanges(seqnames=new_low_ref[,1], ranges=IRanges::IRanges(start=new_low_ref[,2], end=new_low_ref[,3]), mcols= new_low_ref[,-c(1:3)])			
			
			het <- c(het,new_low_het)
			hom <- c(hom, new_low_hom)
                }

		if (adjust_method != "minimal") {

			#Merging inversions and checking whether they have the same genotype
			new_regions <- list(GenomicRanges::reduce(het), GenomicRanges::reduce(hom))
        	       	new_inversions <- genotype_new_regions(new_regions, het, hom, reads, confidence, base, sex, prior, prior_male)	
		
		} else {
	
			#Trying to find and genotype the intersection of all inversions in a cluster of overlaps
			new_regions <- list(minimal(het), minimal(hom))
                        new_inversions <- genotype_new_regions(new_regions, het, hom, reads, confidence, base, sex, prior, prior_male)
			
		}



		if (adjust_method == "deltas" | adjust_method == "all" | adjust_method == "low") {
			
			#If there are still overlapping inversions after trying adjustment and merging, just take the largest.
			het <- GenomicRanges::sort(new_inversions[[1]], by = ~width, decreasing=TRUE)
			hom <- GenomicRanges::sort(new_inversions[[2]], by = ~width, decreasing=TRUE)
			het_merged <- GenomicRanges::reduce(het)
			hom_merged <- GenomicRanges::reduce(hom)
			het <- het[GenomicRanges::findOverlaps(het_merged, het, select="first")]
                        hom <- hom[GenomicRanges::findOverlaps(hom_merged, hom, select="first")]
		
			#Combining the inversions from the different adjustment methods
			all <- c(het,hom)
		
		} else {

			#For "merge" and "minimal" don't bother about the ones that couldn't be merged
	                all <- c(new_inversions[[1]], new_inversions[[2]])

		}

		if (adjust_method == "low" | adjust_method == "all"){


			no_inversion <- inversions[inversions$probability<confidence,]
			no_inversion <- GenomicRanges::GRanges(seqnames=no_inversion[,1], ranges=IRanges::IRanges(start=no_inversion[,2], end=no_inversion[,3]), mcols= no_inversion[,-c(1:3)])
			no_inversion <- no_inversion[is.na(GenomicRanges::findOverlaps(no_inversion, all, select="first")),]
			no_inversion <- no_inversion[is.na(GenomicRanges::findOverlaps(no_inversion, new_low_ref, select="first")),]
			no_inversion <- data.frame(chr=GenomicRanges::seqnames(no_inversion), start=GenomicRanges::start(no_inversion), end=GenomicRanges::end(no_inversion), GenomicRanges::mcols(no_inversion))
			new_low_ref <- data.frame(chr=GenomicRanges::seqnames(new_low_ref), start=GenomicRanges::start(new_low_ref), end=GenomicRanges::end(new_low_ref), GenomicRanges::mcols(new_low_ref))

		} else {
			
			no_inversion <- inversions[inversions$probability<confidence,]	

		}


		#Constructing a dataframe of the adjusted inversions and the reference homozygotes
               	all <- data.frame(chr=GenomicRanges::seqnames(all), start=GenomicRanges::start(all), end=GenomicRanges::end(all), GenomicRanges::mcols(all))
               	ref <- inversions[inversions$probability>=confidence &  ( inversions$genotype == "0|0" | inversions$genotype == "0"),]
		colnames(all) <- colnames(ref)
		colnames(no_inversion) <- colnames(ref)
		colnames(new_low_ref) <- colnames(ref)			
                inversions <- rbind(all, no_inversion, ref, new_low_ref)

                return(inversions)

        }


