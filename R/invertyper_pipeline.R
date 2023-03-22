#' The BED files assume 1-based coordinates.
#'
#' @param regions_to_genotype
#' @param prior
#' @param prior_male
#' @param adjust_method
#' @param inputfolder
#' @param outputfolder
#' @param sex
#' @param vcf
#' @param paired_reads
#' @param confidence
#' @param blacklist
#' @param chromosomes
#' @param numCPU
#' @param save_composite_files
#' @param write_browser_files
#' @param discover_breakpointr_inversions
#' @param breakpointr_prior
#' @param breakpointr_prior_male
#' @param windowsize
#' @param minReads
#' @param background
#' @param output_file
#'
#' @return  [...]
#'
#' @export

invertyper_pipeline <- function(regions_to_genotype=NULL, prior = c(0.333,0.333,0.333), prior_male=c(0.5,0.5),  adjust_method=c("raw","merge","deltas","minimal", "low", "all"), inputfolder='./', outputfolder='./', sex='female', vcf=NULL, paired_reads=TRUE, confidence=0.95, blacklist=NULL, chromosomes=NULL, numCPU=24, save_composite_files=FALSE, write_browser_files=FALSE, discover_breakpointr_inversions=TRUE, breakpointr_prior=c(0.9,0.05,0.05), breakpointr_prior_male=c(0.9,0.1), windowsize=c(40,120,260), minReads=c(15,50,50), background=0.2, output_file="inversions.txt") {

stopifnot("The provided sex should be 'male' or 'female'. This affects whether we consider the X chromosome diploid, and whether we look for inversions on a haploid Y chromosome" = (sex %in% c('male', 'female') & length(sex) ==1))
stopifnot("There is nothing to genotype: either provide a list of putative inversions (regions_to_genotype) or set discover_breakpointr_inversions=TRUE" = !is.null(regions_to_genotype) | discover_breakpointr_inversions)
stopifnot("The confidence threshold for posterior probabilties must be between 0 and 1 (and it should be of type numeric). In general, the arbitrary convention is to choose 0.95." = (as.numeric(confidence) < 1 & as.numeric(confidence) > 0 & is.numeric(confidence)))
stopifnot('Please choose a value for adjust_method. This controls how InvertypeR will attempt to improve the inversion coordinates you supply. Valid values are "raw","merge","deltas","minimal", "low", or "all".' = all(length(adjust_method)==1 & adjust_method %in% c("raw","merge","deltas","minimal", "low", "all")) | is.null(regions_to_genotype))
stopifnot("The arguments minReads and windowsize should be parallel (have the same length)" = (length(minReads)==length(windowsize)))

if(all(prior == c(0.333,0.333,0.333)) | (sex=='male' & all(prior_male == c(0.5,0.5)))){

message("Using the default priors (prior for females, or both prior and prior_male for males) is not recommended. Consider what fraction of the putative inversions you wish to genotype are likely to have non-reference genotypes.")
warning("Using the default priors (prior for females, or both prior and prior_male for males) is not recommended. Consider what fraction of the putative inversions you wish to genotype are likely to have non-reference genotypes.") 

}



if(!is.null(regions_to_genotype)){
	regions_to_genotype <- import(regions_to_genotype)
}

if(!is.null(blacklist)){
	blacklist <- import(blacklist)
}

composite_files <- create_composite_files(inputfolder=inputfolder, type=c("wc","ww"), numCPU=numCPU, vcf=vcf, paired_reads=paired_reads, blacklist=blacklist, outputfolder=outputfolder, save_composite_files=save_composite_files, 
chromosomes=chromosomes)
composite_files[[1]] <- drop_seq_qual(composite_files[[1]], paired_reads=paired_reads)
composite_files[[2]] <- drop_seq_qual(composite_files[[2]], paired_reads=paired_reads)

if(discover_breakpointr_inversions){
message('start inv discovery')
possible_inversions <- discover_possible_inversions(composite_files=composite_files, windowsize=windowsize, minReads=minReads, paired_reads=paired_reads, numCPU=numCPU, chromosomes=chromosomes, blacklist=blacklist, background=background)

message('start inv genotyping (bpr)')
breakpointr_inversions <- invertyper(WW_reads=composite_files$WW, WC_reads=composite_files$WC, regions_to_genotype=possible_inversions, blacklist=blacklist, paired_reads=paired_reads, sex=sex, confidence=confidence,prior=breakpointr_prior, prior_male=breakpointr_prior_male, output_file=paste0("breakpointr_",output_file), adjust_method=c("merge"))




if( write_browser_files){

write_UCSC_browser_files(inversions=breakpointr_inversions, WW_reads=composite_files$WW, WC_reads=composite_files$WC, confidence=confidence, paired_reads=paired_reads, outputfolder=outputfolder, prefix='breakpointr_')

}


}
message('now for std genotyping')



if(!is.null(regions_to_genotype)){

inversions <- invertyper(WW_reads=composite_files$WW, WC_reads=composite_files$WC, regions_to_genotype=regions_to_genotype, blacklist=blacklist, paired_reads=paired_reads, sex=sex, confidence=confidence,prior=prior,prior_male=prior_male, output_file=output_file, adjust_method=adjust_method)



if( write_browser_files){

write_UCSC_browser_files(inversions=inversions, WW_reads=composite_files$WW, WC_reads=composite_files$WC, confidence=confidence, paired_reads=paired_reads, outputfolder=outputfolder, prefix='')

}

}

message('about to return for the last time')





if(!is.null(regions_to_genotype) & discover_breakpointr_inversions){

return(list(inversions, breakpointr_inversions))

} else if(!is.null(regions_to_genotype)) {

return(inversions)

} else {

return(breakpointr_inversions)

}

}
