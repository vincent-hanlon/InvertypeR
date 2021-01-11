#' Counts reads in genomic regions
#'
#' For each composite file, this function tallies up the forward and reverse reads in each putative inversion (region).
#'
#' One of the slowest steps of the genotyping process when composite files are large.
#'
#' @param WW_bam A GRanges object for the WW composite file
#' @param WC_bam A GRanges object for the WC composite file
#' @param regions A GRanges object for putative inversions
#' @return A dataframe with read counts by file
#'
#'
#' @export
count_regions <- function(WW_reads, WC_reads, regions) {

       #Extracting counts of W and C reads by region
        WW_w <- GenomicRanges::countOverlaps(subject=WW_reads[GenomicRanges::strand(WW_reads)=="+"],query=regions)
        WW_c <- GenomicRanges::countOverlaps(subject=WW_reads[GenomicRanges::strand(WW_reads)=="-"],query=regions)
        WC_w <- GenomicRanges::countOverlaps(subject=WC_reads[GenomicRanges::strand(WC_reads)=="+"],query=regions)
        WC_c <- GenomicRanges::countOverlaps(subject=WC_reads[GenomicRanges::strand(WC_reads)=="-"],query=regions)

        #Assemmbling a data frame with region info
        counts <- data.frame(GenomicRanges::seqnames(regions), GenomicRanges::start(regions), GenomicRanges::end(regions), WW_w, WW_c, WC_w, WC_c)
        colnames(counts) <- c("chr","start","end", "WW_w", "WW_c", "WC_w","WC_c")

        return(counts)
}






#' Reads BAM files into R
#'
#' Reads near or within putative inversions (regions) that do not overlap the blacklist are read into GRanges objects.
#'
#' Requires mapQ>=10 and excludes unmapped and duplicate reads.
#'
#' @param WW_bam Path to the WW composite files
#' @param WC_bam Path to the WC composite file
#' @param regions A GRanges object for putative inversions
#' @param paired_reads A boolean: are the reads paired-end or not? Default TRUE.
#' @param blacklist Path to a BED file of regions to exclude
#' @return A list containing two GRanges objects of reads
#'
#'
#' @export
read_regions <- function(WW_bam, WC_bam, regions, paired_reads=TRUE, blacklist=blacklist) {

	if(blacklist!=""){
        	regions <- GenomicRanges::setdiff(GenomicRanges::reduce((regions+1000000),min.gapwidth=1000),import(blacklist))
        } else {
        	regions <- GenomicRanges::reduce((regions+1000000),min.gapwidth=1000)
        }


        if(paired_reads) {

		#Reading in files. Warnings suppressed because GenomicAlignments::readGAlignment... generates them for no good reason
		#Note that the intervals are extended by 1Mb (so that we can read in flanking regions for other purposes) and reduced so they aren't overlapping
                p1 <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired=TRUE,isUnmappedQuery=FALSE,isDuplicate=FALSE),
			mapqFilter=10,which=regions,what=c("mapq"))
                WW_reads <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(WW_bam,param=p1))
                WC_reads <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(WC_bam,param=p1))

        } else {

                p1 <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired=FALSE,isUnmappedQuery=FALSE,isDuplicate=FALSE),mapqFilter=10,which=regions,what=c("mapq"))
                WW_reads <- suppressWarnings(GenomicAlignments::readGAlignments(WW_bam,param=p1))
                WC_reads <- suppressWarnings(GenomicAlignments::readGAlignments(WC_bam,param=p1))

         }

return(list(WW_reads,WC_reads))

}
 




#' Mandatory preprocessing of reads for BreakpointR (Porubsky et al. 2020)
#'
#' Code from BreakpointR, https://doi.org/10.1093/bioinformatics/btz681
#'
#' To use the deltaWCalculator from BreakpointR with paired-end reads, they must first be merged into one fragement. It's not supplied as a function so we have to reproduce a piece of the code.
#'
#' @param data.raw A GAlignmentPairs object
#' @return A GenomicRanges object
#'
#'
#' @export
pair2frgm <- function(data.raw){

        data.prop.pairs <- data.raw[GenomicAlignments::isProperPair(data.raw)]
        data.first <- as(GenomicAlignments::first(data.prop.pairs), 'GRanges')
        data.last <- as(GenomicAlignments::last(data.prop.pairs), 'GRanges')

        data.first.plus <- data.first[GenomicRanges::strand(data.first) == '+']
        data.first.minus <- data.first[GenomicRanges::strand(data.first) == '-']
        data.last.plus <- data.last[GenomicRanges::strand(data.last) == '+']
        data.last.minus <- data.last[GenomicRanges::strand(data.last) == '-']

        frag.plus.mapq <- data.first.plus$mapq + data.last.minus$mapq
        frag.minus.mapq <- data.first.minus$mapq + data.last.plus$mapq

        data.frag.plus <- GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(data.first.plus), ranges=IRanges::IRanges(start=GenomicRanges::start(data.first.plus), end=GenomicRanges::end(data.last.minus)), strand=GenomicRanges::strand(data.first.plus), mapq=frag.plus.mapq)
        GenomeInfoDb::seqlengths(data.frag.plus) <- GenomeInfoDb::seqlengths(data.first)
        data.frag.minus <- GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(data.first.minus), ranges=IRanges::IRanges(start=GenomicRanges::start(data.last.plus), end=GenomicRanges::end(data.first.minus)), strand=GenomicRanges::strand(data.first.minus), mapq=frag.minus.mapq)
        GenomeInfoDb::seqlengths(data.frag.minus) <- GenomeInfoDb::seqlengths(data.first)


        data <- GenomicRanges::sort(c(data.frag.plus, data.frag.minus), ignore.strand=TRUE)


        return(data)
}




#' Read processing for inversion adjustment
#'
#' Input reads are sent through the deltaWCalculator from BreakpointR, which associates them with deltaW values approximating the local change in read direction  
#'
#' For each read, we sum multiple deltaW values calculated for binsizes 5, 10, 10, 20, 40, and 80. Yes, 10 is there twice: it means we're prioritizing changes in read direction from
#' the smaller scale. 
#'
#' @param WW_reads A GRanges object from the WW composite file.
#' @param WC_reads A GRanges object from the WC composite file.
#' @param paired_reads Boolean: are the reads paired-end? Default TRUE.
#' @return A list of four GRanges objects: the input objects WW_reads and WC_reads and the two output objects (for WW and WC composite files) that have been annotated with deltaW values.
#'
#'
#' @export
process_reads <- function(WW_reads, WC_reads, paired_reads=TRUE){

        #Processing input reads for Aaron's deltaWCalculator

        if (paired_reads) {

                WW_reads <- pair2frgm(WW_reads)
                WC_reads <- pair2frgm(WC_reads)


        } else {

                WW_reads <- GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(WW_reads),ranges=IRanges::IRanges(start=GenomicRanges::start(WW_reads), end=GenomicRanges::end(WW_reads)), strand=GenomicRanges::strand(WW_reads), GenomicRanges::mcols(WW_reads), seqlengths=GenomeInfoDb::seqlengths(WW_reads))
                WC_reads <- GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(WC_reads),ranges=IRanges::IRanges(start=GenomicRanges::start(WC_reads), end=GenomicRanges::end(WC_reads)), strand=GenomicRanges::strand(WC_reads), GenomicRanges::mcols(WC_reads), seqlengths=GenomeInfoDb::seqlengths(WC_reads))

        }

        #Calculating deltaW values
        WW_d <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 5))
        WC_d <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 5))

        WW_d_10 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 10))
        WC_d_10 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 10))

        WW_d_20 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 20))
        WC_d_20 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 20))

        WW_d_40 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 40))
        WC_d_40 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 40))

        WW_d_80 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 80))
        WC_d_80 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 80))

        #Taking an arbitrary sum of deltaWs to look near and far (including 2x the deltaW10)
        WW_d$deltaW <- WW_d$deltaW + 2*WW_d_10$deltaW + WW_d_20$deltaW + WW_d_40$deltaW + WW_d_80$deltaW
        WC_d$deltaW <- WC_d$deltaW + 2*WC_d_10$deltaW + WC_d_20$deltaW + WC_d_40$deltaW + WC_d_80$deltaW



        return(list(WW_reads, WC_reads,WW_d,WC_d))

}


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
gr <- GenomicRanges::GRanges(seqnames=bed[,1], ranges=IRanges::IRanges(start=bed[,2], end=bed[,3]))
return(gr)

}

