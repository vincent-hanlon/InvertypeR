
#' Checks if adjusted inversions have the same genotype
#'
#' Runs genotype_inversions(), discards inversions that don't have the same genotype as before adjustment, and then adds back in the originals for the ones that weren't adjusted.
#'
#' @param new_regions A list of length two: adjust_deltaW() output for the het and hom cases.
#' @param het A GRanges object of heterozygous inversions.
#' @param hom A GRanges object of homozygous inversions.
#' @param reads reads A list of two Granges objects: reads from a WW composite file and a WC composite file.
#' @param confidence Posterior probability threshold above which you consider genotype calls to be reliable. Used to decide whether to keep adjusted inversions. Default 0.95.
#' @param base A list output by WWCC_background().
#' @param sex Sex of sample to figure out sex chromosomes. Default "female".
#' @param prior Vector of three prior weights for inversion genotypes. For example, c("ref","het","hom") = c(0.9,0.05,0.05).
#' @param prior_male Vector of two prior weights for male sex chromosomes. For example, c("ref", "inv") = c(0.9,0.1). 
#' @param paired_reads Boolean: are the reads paired-end? Default TRUE.
#' @return A list containing two GRanges objects: the adjusted heterozygous inversions and the adjusted homozygous inversions.
#'
#'
#' @export
genotype_new_regions <- function(new_regions, het, hom, reads, confidence, base,  sex="female", prior, prior_male, paired_reads=TRUE){

		#Re-genotyping the adjusted inversions
               	new_het <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=new_regions[[1]], background=base[[1]], base_state=base[[2]],  sex=sex,
                        prior=prior, prior_male=prior_male)
               	new_hom <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=new_regions[[2]], background=base[[1]], base_state=base[[2]],  sex=sex,
                        prior=prior, prior_male=prior_male)


		#Discarding the ones that don't match genotypes and constructing a dataframe
		new_het <- new_het[(new_het$genotype=="0|1" | new_het$genotype =="1|0") & new_het$probability >= confidence,]
		new_hom <- new_hom[(new_hom$genotype=="1|1" | new_hom$genotype =="1") & new_hom$probability >= confidence,]
		new_het <- GenomicRanges::GRanges(seqnames=new_het[,1], ranges=IRanges::IRanges(start=new_het[,2], end=new_het[,3]), mcols= new_het[,-c(1:3)])
		new_hom <- GenomicRanges::GRanges(seqnames=new_hom[,1], ranges=IRanges::IRanges(start=new_hom[,2], end=new_hom[,3]), mcols=new_hom[,-c(1:3)])

		#Retrieving inversions that weren't successfully adjusted
		het <- het[is.na(GenomicRanges::findOverlaps(het, new_het, select="first")),]
                hom <- hom[is.na(GenomicRanges::findOverlaps(hom, new_hom, select="first")),]

		return(list(c(het,new_het), c(hom,new_hom)))

}
