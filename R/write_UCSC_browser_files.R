#' Writes genome browser files for inversions
#'
#' Given a list of inversions from invertyper() and two composite files, it makes nicely-formatted BED files compatible with the UCSC Genome Browser
#' The inversions are coloured by genotype (black is 1/1, grey is 0/1, light grey is 0/0), and the Strand-seq reads are coloured by strand as in BreakpointR etc.
#'
#' @param inversions
#' @param WW_reads
#' @param WC_reads
#' @param confidence
#' @param paired_reads
#' @param outputfolder
#' @param prefix
#'
#' @return  [...]
#'
#' @export

write_UCSC_browser_files <- function(inversions=NULL, WW_reads=NULL, WC_reads=NULL, confidence=0.95, paired_reads=TRUE, outputfolder="./", prefix=''){

if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
}


region <- inversions[inversions[,9]>confidence & inversions[,8]!=0 & inversions[,8]!="0|0",]

if(nrow(region)>0){

region <- GenomicRanges::reduce(widen(granges=GenomicRanges::makeGRangesFromDataFrame(region, ignore.strand=TRUE, seqnames.field=names(region)[1], start.field=names(region)[2],end.field=names(region)[3]), seqlengths=GenomeInfoDb::seqlengths(WW_reads), distance=5e06), min.gapwidth=1000)

WW_reads <- galignment_to_granges(WW_reads, paired_reads=paired_reads, purpose='BreakpointR', region=region)
WC_reads <- galignment_to_granges(WC_reads, paired_reads=paired_reads, purpose='BreakpointR', region=region)

suppressMessages(breakpointR::breakpointr2UCSC(paste0(prefix,"WW_composite_file"), outputDirectory=outputfolder, fragments=WW_reads))
suppressMessages(breakpointR::breakpointr2UCSC(paste0(prefix,"WC_composite_file"), outputDirectory=outputfolder, fragments=WC_reads))

inversions <- inversions[inversions[,9]>confidence,]
inversions[,2] <- inversions[,2]-1
inversions[,3] <- inversions[,3]-1 
inversions[,4] <- inversions[,8]
inversions[,5] <- 0
inversions[,6] <- "."
inversions[,7] <- inversions[,2]
inversions[,8] <- inversions[,2]
inversions[inversions[,4]==0 | inversions[,4]=="0|0",9] <- "250,250,250"
inversions[inversions[,4]=="1|0" | inversions[,4]=="0|1",9] <- "125,125,125" 
inversions[inversions[,4]=="1|1" | inversions[,4]=="1",9] <- "0,0,0"
inversions[,10] <- NULL

savefile.invs <- file.path(outputfolder, paste0(prefix,'inversions.bed.gz'))
savefile.invs.gz <- gzfile(savefile.invs, 'w') 
utils::write.table(paste0("track name=",prefix,"inversions itemRgb=On"), file=savefile.invs.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE, sep='\t')
utils::write.table(inversions,file=savefile.invs.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
close(savefile.invs.gz)
} else {

 warning("There are no confident inversions to create a UCSC browser file")

}


}
