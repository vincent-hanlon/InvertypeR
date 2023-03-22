
#' Creates composite files for inversion genotyping
#'
#' 
#'
#' @param inputfolder
#' @param type
#' @param numCPU
#' @param vcf
#' @param paired_reads
#' @param blacklist
#' @param outputfolder
#' @param save_composite_files
#' @param chromosomes
#'
#' @return  [...]
#'
#' @export


#' This sometimes has issues with memory usage: if the memory usage on the first BPR step gets too high, R will crash.
#' be sure to include in the documentation that the export function from the rtracklayer package can be used to write a GenomicAlignment or GenomicAlignmentPairs composite file to a BAM file, which may be useful for bioinformatic command line tools.
create_composite_files <- function(inputfolder='./', type=c("wc","ww"), numCPU=24, vcf=NULL, paired_reads=TRUE, blacklist=NULL, outputfolder='./', save_composite_files=FALSE, chromosomes=NULL){

message('start composite')

stopifnot("The type argument should be 'wc' for a Watson-Crick composite file, 'ww' for a Watson-Watson composite file, or c('ww','wc') for both." = (!any(!type %in% c('wc','ww'))))
stopifnot("The path to a VCF file is required for a Watson-Crick ('wc') composite file." = !(is.null(vcf) & 'wc' %in% type))

if(is.null(blacklist)){
warning("A blacklist is highly recommended both for creating composite files and for genotyping inversions.")
}

if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
}

bamlist <- list.files(inputfolder, pattern="\\.bam$")


suppressMessages(cl <- parallel::makeCluster(numCPU))
#(maybe should include something here like this:  doParallel::registerDoParallel(cl))
# We first must read in the BAM files without removing reads in blacklisted regions or chromosomes that are left out
message('start read bams')
galignmentslist <- parallel::parSapply(cl, bamlist, read_bam, blacklist=NULL, chromosomes=NULL, paired_reads=paired_reads)
names(galignmentslist) <- bamlist
message('start granges conversion')
grangeslist <- parallel::parSapply(cl, galignmentslist, galignment_to_granges, purpose='BreakpointR', paired_reads=paired_reads, pair2frgm=FALSE, simplify = FALSE,USE.NAMES = TRUE)
names(grangeslist) <- names(galignmentslist)
parallel::stopCluster(cl)
message('start bpr')
# breakpointR unfortunately can't find the strand state of segments in the genome if reads in blacklisted regions are first removed. 
# This is probably because there aren't enough Strand-switches
# Could be a github issue: subtract blacklisted regions from BAM file directly?
# This means it's best to read in the offending reads, and then only remove them later on.
# ALSO: the galignments_to_granges should really just take the first read for PE reads, unless pair2frgm is specified
bpr <- suppressMessages(breakpointr_for_invertyper(grangeslist, paired_reads=FALSE, numCPU=8, windowsize=20000000, binMethod="size", minReads=50, background=0.2, maskRegions=blacklist, chromosomes=NULL))

message('done bpr')

rm('grangeslist')
invisible(gc())

message('start chromosome downsample')
if(!is.null(blacklist) | !is.null(chromosomes)){
suppressMessages(cl <- parallel::makeCluster(numCPU))

# now we "re-read" the BAM files, but in reality we're just subsetting Galignments objects according to the chromosomes and blacklist
galignmentslist <- parallel::parSapply(cl, galignmentslist, read_bam, blacklist=blacklist, chromosomes=chromosomes)
names(galignmentslist) <- bamlist

parallel::stopCluster(cl)
}

message('start extract regions')

WWCCregions <- sort(find_regions_with_strand_state(bpr, states=c('ww','cc'), region_size=100000))
WWCCregions$Cs <- NULL
WWCCregions$Ws <- NULL
WWCCregions$state <- NULL
GenomeInfoDb::seqlevels(WWCCregions) <- GenomeInfoDb::seqlevels(galignmentslist[[1]])

if('wc' %in% type){

WCregions <- sort(find_regions_with_strand_state(bpr, states=c('wc'), region_size=100000))
WCregions$Cs <- NULL
WCregions$Ws <- NULL
WCregions$state <- NULL
GenomeInfoDb::seqlevels(WCregions) <- GenomeInfoDb::seqlevels(galignmentslist[[1]])
rm('bpr')
invisible(gc())
message('start phasing')
all_phased_WCregions <- strandPhaseR_for_invertyper(numCPU=numCPU, positions=vcf, WCregions=WCregions, chromosomes=chromosomes, paired_reads=paired_reads, num.iterations=3, galignmentslist=galignmentslist)

} else{
rm('bpr')
invisible(gc())
}

message('start composite assembly')




# maybe this should eventually be parallelized using R parallel and parSapply, but to be honest parallel::mclapply seems faster and less complicated so far.
#cl <- parallel::makeCluster(numCPU)
# parallelizing hides error messages in sub-functions though: so these should be in try-catch statements with useful error messages
# ... and definitely the numCPU=1 case should return useful messages (for example, maybe it needs to be wrapped in an if statement where the alternative is a for loop, like in DP's code.

CC <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state, galignmentslist=galignmentslist, regions=WWCCregions, paired_reads=paired_reads, states='cc', flip_reads=TRUE, mc.cores=numCPU)
WW <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state, galignmentslist=galignmentslist, regions=WWCCregions, paired_reads=paired_reads, states='ww', flip_reads=FALSE, mc.cores=numCPU)

if('wc' %in% type){
message('choose wc composite reads')
CW <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state, galignmentslist=galignmentslist, regions=all_phased_WCregions, paired_reads=paired_reads, states='cw', flip_reads=TRUE, mc.cores=numCPU)
WC <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state, galignmentslist=galignmentslist, regions=all_phased_WCregions, paired_reads=paired_reads, states='wc', flip_reads=FALSE, mc.cores=numCPU)

}

rm('galignmentslist')
invisible(gc())


#parallel::stopCluster(cl)
message('make ww file')
WW_composite_file <- do.call('c',c(CC,WW))

message('make ww file to plot')
to_plot <- list(galignment_to_granges(WW_composite_file[sort(sample(length(WW_composite_file),min(length(WW_composite_file),500000)))], purpose='BreakpointR', paired_reads=paired_reads,  pair2frgm=FALSE))
names(to_plot)[1] <- 'WW_composite_file'

if('wc' %in% type){
message('assemble wc file')

WC_composite_file <- do.call('c',c(WC,CW))
message('make wc toplot')
to_plot <- append(to_plot,galignment_to_granges(WC_composite_file[sort(sample(length(WC_composite_file),min(length(WC_composite_file),500000)))], purpose='BreakpointR', pair2frgm=FALSE, paired_reads=paired_reads))

names(to_plot)[2] <- 'WC_composite_file'

}







message('start plot composites')
invisible(suppressMessages(breakpointr_for_invertyper(to_plot, paired_reads=FALSE, plotspath=outputfolder,numCPU=numCPU, windowsize=1000000, binMethod="size", minReads=50, background=0.2,maskRegions=blacklist,chromosomes=chromosomes)))

if(save_composite_files){

	if('ww' %in% type){
		save(WW_composite_file,file=file.path(outputfolder,"WW_composite_file.RData"))
	}

	if('wc' %in% type){
		save(WC_composite_file,file=file.path(outputfolder,"WC_composite_file.RData"))
	}
}



composite <- list(WW_composite_file, WC_composite_file)
names(composite) <- c("WW", "WC") 


return(composite)

}


