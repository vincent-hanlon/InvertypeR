#' Main wrapper function for the InvertypeR package.
#'
#' Genotypes inversions and adjusts the start and end coordinates as specified by the user. 
#'
#' Prepare the composite files using the scripts available at https://github.com/vincent-hanlon/InvertypeR. The BED files assume 1-based coordinates. We highly recommend using a blacklist,
#' and we provide one on the GitHub page. This package has only been tested for paired-end reads so far. Priors should chosen to reflect the number of inversions you expect to find in an individual (up to
#' ~100 if you use, for example, all publised inversion coordinates (2020)) and the number of unique intervals in your BED file of putative inversions (many ways to calculate this, but 80% overlapping 
#' intervals could be considered the same, for example). If you are confident in the coordinates of your putative inversions and you're just starting out, use adjust_method="raw". Otherwise, use 
#' adjust_method="deltas" or "all". The other options are specific to particular use cases. 
#'
#' Note that we can phase inversions only because the WC composite file already has phased reads. This means that we know a 0|1 inversion on chr1 is on the same homolog as all 0|1 inversions
#' on chr1 in the sample, and that all chr1 1|0 inversions are on the other homolog. However, we don't know whether a 0|1 inversion on chr1 and a 0|1 inversion chr2 came from the same parent.
#' 0|1 inversions are distinguished from 1|0 inversions based on the strand switch in the WC composite file ( WC -> WW or WC -> CC).
#'
#' @param WW_bam Path to a WW composite BAM file containing Strand-seq data.
#' @param WC_bam Path to a WC composite BAM file containing Strand-seq data.
#' @param bed Path to a BED file containing genomic intervals to be genotyped (putative inversions).
#' @param blacklist Path to a BED file containing regions genomic intervals with poor-quality Strand-seq data. Reads that overlap these intervals will not be used. Highly recommended. Default "".
#' @param paired_reads Boolean: are the reads paired-end? Single-end reads have not yet been tested. Default TRUE. 
#' @param sex Sex of sample to figure out sex chromosomes. Default "female".
#' @param confidence Posterior probability threshold above which you consider genotype calls to be reliable. Used to decide whether to keep adjusted inversions, as well as to identify low-confidence 
#'   calls for adjust_method "low". Default 0.95.
#' @param prior Vector of three prior weights for inversion genotypes. For example, c("REF","HET","HOM") = c(0.9,0.05,0.05). Default c(0.33,0.33,0.33).
#' @param prior_male Vector of two prior weights for inversions on the male sex chromosomes. For example, c("REF", "INV") = c(0.9,0.1). Default c(0.5,0.5).
#' @param output_file Name of the file to write to. Default "inversions.txt".
#' @param adjust_method One of "raw", "merge", "deltas", "minimal", "low", or "all". Specifies which method to use to adjust the inversion coordinates (start- and end-points). The adjustment routine 
#'   ensures that adjusted inversions have the same genotype as they did before adjustment (if applicable), and that they overlap at least one of the original unadjusted inversions in every cluster of 
#'   overlapping events. If after 
#'   adjustment (except with "raw") overlapping inversions still remain, we merge the ones of the same genotype. If there are still overlapping inversions, we take the largest.
#'   "raw": No adjustment. "merge": Merge overlapping confident (i.e., posterior probability >= confidence for a genotype) inversions of the same genotype. "Deltas": adjust the endpoints of confident 
#'   inversions based on the Strand-seq data. 
#'   We calculate deltaW values for each read, which measure the change in read direction for nearby reads, and take the coordinates that correspond to the two highest values (i.e. the two spots with 
#'   the greates change in read direction, which we hope will correspond to inversion breakpoints). This is done once for every 
#'   set of overlapping inversions of the same genotype. "minimal": An experimental method. For confident overlapping inversions of 
#'   the same genotype, take the interval common to all of them if one exists. "low": Much like "deltas", except that we only adjust inversions for which no confident genotype could be found (i.e., all 
#'   prior probabilities < confidence). Sometimes this will "find" a confident inversion if the original coordinates were slightly off. "all": Does both "deltas" and "low".
#' @return A dataframe that associates each entry in the BED file (putative inversions) with a genotype and posterior probability. Cols 1-3 are genomic coordinates, and cols 4-7 are read counts for the 
#'   two composite files (WC and WW, or WC and CC if that's what you submitted) subset by Watson (W) and Crick (C) strands. Col 8 is the most probable genotype, and col 9 is the posterior probability 
#'   associated with that genotype. Col 10 labels inversions with low read density. The start and end coordinates of such inversions should be checked by viewing Strand-seq data in a genome browser.
#'
#' @export
invertyper <- function(WW_bam, WC_bam, bed, blacklist="", paired_reads=TRUE, sex="female", confidence=0.95,prior=c(0.333,0.333,0.333), prior_male=c(0.5,0.5), 
					output_file="inversions.txt", adjust_method=c("raw","merge","deltas","minimal", "low", "all")){

	#Loading the BED files as GRanges objects
	regions <- import(bed)

	#Takes reads from the BAM files that overlap the regions and stores them as GRanges
	reads <- read_regions(WW_bam, WC_bam, regions, paired_reads=paired_reads, blacklist=blacklist)

        #Accurate background estimate, plus base strand state for the WW/CC file
        base <- WWCC_background(WW_bam, binsize=1000000, paired_reads=paired_reads)

	#Genotyping the inversions
	inversions <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=regions, background=base[[1]], base_state=base[[2]],  sex=sex, 
		prior=prior, prior_male=prior_male)

	#two methods under development for dealing with overlapping intervals
	if( adjust_method != "raw" ){
		
		inversions <- adjust_intervals(inversions=inversions, reads=reads, confidence=confidence, base=base,  sex=sex, prior=prior, prior_male=prior_male, 
			adjust_method=adjust_method, paired_reads=paired_reads)
	}

	#Marking inversions with <0.001 reads per bp as being low-density, because they may span centromeres or chromosome ends/
	#The 0.001 figure is somewhat arbitrary, but it's adjusted based on the read count in the WW file of any submitted inversions. 
	inversions[inversions[,9]>=confidence & inversions[,8]!=0 & inversions[,8]!="0|0" & rowSums(inversions[,c(4:7)])/(inversions[,3]-inversions[,2] + 1) < (0.001*base[[3]]/27734969), 10] <- 
			"low read density: check for inaccurate inversion coordinates"

	colnames(inversions)[10] <- "low_read_density"

	write.table(inversions, output_file, quote=FALSE, sep="\t", row.names=FALSE)

	return(inversions)
}

