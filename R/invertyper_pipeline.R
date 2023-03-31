#' The BED files assume 1-based coordinates.
#'
#' @param regions_to_genotype
#' @param prior
#' @param haploid_prior
#' @param adjust_method
#' @param inputfolder
#' @param outputfolder
#' @param haploid_chromosomes
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
#' @param breakpointr_haploid_prior
#' @param windowsize
#' @param minReads
#' @param background
#' @param output_file
#'
#' @return  [...]
#'
#' @export

invertyper_pipeline <- function(regions_to_genotype=NULL, prior = c(0.333,0.333,0.333), haploid_prior=c(0.5,0.5),  adjust_method=c("raw","merge","deltas","minimal", "low", 
"all"), inputfolder='./', outputfolder='./', haploid_chromosomes=NULL, vcf=NULL, paired_reads=TRUE, confidence=0.95, blacklist=NULL, chromosomes=NULL, numCPU=24, 
save_composite_files=FALSE, write_browser_files=FALSE, discover_breakpointr_inversions=TRUE, breakpointr_prior=c(0.9,0.05,0.05), breakpointr_haploid_prior=c(0.9,0.1), windowsize=c(40,120,360), minReads=c(15,50,50), background=0.2, output_file="inversions.txt") {

stopifnot("There is nothing to genotype: either provide a list of putative inversions (regions_to_genotype) or set discover_breakpointr_inversions=TRUE" = !is.null(regions_to_genotype) | discover_breakpointr_inversions)
stopifnot("The confidence threshold for posterior probabilties must be between 0 and 1 (and it should be of type numeric). In general, the arbitrary convention is to choose 0.95." = (as.numeric(confidence) < 1 & as.numeric(confidence) > 0 & is.numeric(confidence)))
stopifnot('Please choose a value for adjust_method. This controls how InvertypeR will attempt to improve the inversion coordinates you supply. Valid values are "raw","merge","deltas","minimal", "low", or "all".' = all(length(adjust_method)==1 & adjust_method %in% c("raw","merge","deltas","minimal", "low", "all")) | is.null(regions_to_genotype))
stopifnot("The arguments minReads and windowsize should be parallel (have the same length)" = (length(minReads)==length(windowsize)))
stopifnot("The regions_to_genotype file cannot be found" = is.null(regions_to_genotype) || file.exists(regions_to_genotype))
stopifnot("The blacklist file cannot be found" = is.null(blacklist) || file.exists(blacklist))
stopifnot("The vcf file cannot be found" = is.null(vcf) || file.exists(vcf))
stopifnot("The inputfolder cannot be found" = file.exists(inputfolder))
stopifnot("Any haploid chromosomes you specify should be include as chromosomes too. This means including chrX and chrY as chromosomes AND as haploid_chromosomes for human males" = all(haploid_chromosomes %in% chromosomes) || is.null(chromosomes))
if(!is.null(regions_to_genotype) && (all(prior == c(0.333,0.333,0.333)) | (!is.null(haploid_chromosomes) & all(haploid_prior == c(0.5,0.5))))){

warning("Using the default priors (prior for homogametic diploids (e.g., human females), haploid_prior for haploids, or both for heterogametic dipoids (e.g., human males))
is not recommended. Consider what fraction of the putative inversions you wish to genotype are likely to have non-reference genotypes.")

}


if(!is.null(regions_to_genotype)){
	regions_to_genotype <- import_bed(regions_to_genotype)
}

if(!is.null(blacklist)){
	blacklist <- import_bed(blacklist)
}


if(chromosomes==haploid_chromosomes){

type <- 'ww'

} else {

type <- c("wc","ww")

}


composite_files <- create_composite_files(inputfolder=inputfolder, type=type, numCPU=numCPU, vcf=vcf, paired_reads=paired_reads, blacklist=blacklist, outputfolder=outputfolder, save_composite_files=save_composite_files, chromosomes=chromosomes)
composite_files[[1]] <- drop_seq_qual(composite_files[[1]], paired_reads=paired_reads)
composite_files[[2]] <- drop_seq_qual(composite_files[[2]], paired_reads=paired_reads)

if(discover_breakpointr_inversions){
message('start inv discovery')
possible_inversions <- discover_possible_inversions(composite_files=composite_files, windowsize=windowsize, minReads=minReads, paired_reads=paired_reads, numCPU=numCPU, chromosomes=chromosomes, blacklist=blacklist, background=background)

message('start inv genotyping (bpr)')
save(possible_inversions, file="possible.RData")
breakpointr_inversions <- invertyper(WW_reads=composite_files$WW, WC_reads=composite_files$WC, regions_to_genotype=possible_inversions, blacklist=blacklist,paired_reads=paired_reads, haploid_chromosomes=haploid_chromosomes, confidence=confidence,prior=breakpointr_prior, haploid_prior=breakpointr_haploid_prior, output_file=paste0("breakpointr_",output_file), adjust_method=c("merge"))




if( write_browser_files){

write_UCSC_browser_files(inversions=breakpointr_inversions, WW_reads=composite_files$WW, WC_reads=composite_files$WC, confidence=confidence, paired_reads=paired_reads, outputfolder=outputfolder, prefix='breakpointr_')

}


}
message('now for std genotyping')



if(!is.null(regions_to_genotype)){

inversions <- invertyper(WW_reads=composite_files$WW, WC_reads=composite_files$WC, regions_to_genotype=regions_to_genotype, blacklist=blacklist, paired_reads=paired_reads, haploid_chromosomes=haploid_chromosomes, confidence=confidence,prior=prior,haploid_prior=haploid_prior, output_file=output_file, adjust_method=adjust_method)



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
