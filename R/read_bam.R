
#' For composite BAM writing, maybe something related to fwrite could be created? Fastest way to write... write to SAM format and then compress? Or use rtracklayer's function
#' Going to expect sorted, indexed BAM files.
#' Maybe for SNV calling, https://bioconductor.org/packages/release/bioc/manuals/VariantTools/man/VariantTools.pdf should be used? 
#' However, that looks clumsy: it won't take GAlingments format. So it's best to ask for an external VCF, and just suggest in the documentation that this can easily be done with bbtools (give command)
#' If a GenomicAlignments or GRanges object is passed instead of the path to a BAM file, this will just try excluding reads that aren't part of the provided regions. This is helpful because
#' That way InvertypeR can be given either a GAlignments composite file (to throw away useless reads) or a BAM composite file (old-school InvertypeR)

#' @param bam
#' @param region
#' @param paired_reads
#' @param blacklist
#' @param chromosomes
#'
#' @return  [...]
#'
#' @export
read_bam <- function(bam=NULL, region=NULL, paired_reads=TRUE, blacklist=NULL, chromosomes=NULL) {

        if(length(region)==0 & is.character(bam)){
        
	        lengths <- Rsamtools::scanBamHeader(bam)[[1]]$targets
                region <- GenomicRanges::GRanges(names(lengths),ranges=IRanges::IRanges(1,lengths))

        } else if(length(region)>0) {
	
		region <- GenomicRanges::reduce(region, min.gapwidth=1000)

	} else if(length(region)==0 & !is.character(bam) & (!is.null(chromosomes) | !is.null(blacklist))){	

		lengths <- GenomeInfoDb::seqlengths(bam)
		region <- GRanges(seqnames=names(lengths), ranges=IRanges::IRanges(start=c(1:length(lengths)), end=lengths))			

	} else {	
	
		warning("read_bam is returning the input reads unaltered because neither a subsetting region nor the path to a BAM file were supplied.")
		return(bam)	
        }

        if(length(blacklist)>0){
                region <- GenomicRanges::setdiff(region,blacklist)
        }

        if(!is.null(chromosomes)){
		region <- region[as.vector(GenomicRanges::seqnames(region)) %in% chromosomes]
        }



	if(is.character(bam)){

	stopifnot("Set paired_reads=FALSE when your BAM file contains single-end reads, and vice versa" = suppressMessages(Rsamtools::testPairedEndBam(bam))==paired_reads)
	

        if(paired_reads) {

		#Reading in files. Warnings suppressed because GenomicAlignments::readGAlignment... generates them for no good reason
		#Note that the intervals were already extended by 1Mb (so that we can read in flanking region for other purposes) and reduced so they aren't overlapping
                p1 <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired=TRUE,isUnmappedQuery=FALSE,isDuplicate=FALSE, isSecondaryAlignment=FALSE, hasUnmappedMate=FALSE, isProperPair=TRUE),
			mapqFilter=10,which=region,what=c('seq', 'qual','mapq','cigar','flag'))
                reads <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bam,param=p1))

        } else {

                p1 <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired=FALSE,isUnmappedQuery=FALSE,isDuplicate=FALSE, isSecondaryAlignment=FALSE),mapqFilter=10,which=region,what=c('seq', 'qual','mapq','cigar','flag'))
                reads <- suppressWarnings(GenomicAlignments::readGAlignments(bam,param=p1))

         }


	} else if(length(region)>0) {
	
		reads <- bam[S4Vectors::queryHits(GenomicAlignments::findOverlaps(bam, region))]
	} else {
	
		stop("The read_bam function requires either the path to a BAM file, or a GenomicAlignments object full of reads AND a GRanges object containing genomic regions for which the corresponding reads will be selected.")

	}

	return(reads)

}



