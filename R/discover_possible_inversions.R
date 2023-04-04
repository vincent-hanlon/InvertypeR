#'
#' For BreakpointR version 1.16.0, it is currently unwise to use this with a blacklist and few reads/chromosomes, because when no reads intersect the blacklist it removes 
#' all reads. This bug will hopefully be fixed in the next version
#' 
#' @param composite_files A list of 1 or 2 composite files, named WC and WW (or just WW for haploids)
#' @param windowsize A vector of integers. Each integer will be used the set the number of reads in a bin for a BreakpointR run.
#' @param minReads A vector of integers parallel to windowsize, giving the minimum number of reads in a bin for it to be considered/genotyped
#' @param paired_reads Boolean. Paired-end reads?
#' @param numCPU Integer. How many threads to use?
#' @param chromosomes A character vector of chromosome names to examine
#' @param blacklist A GRanges object containing intervals with suspected poor-quality Strand-seq data. 
#' @param background A parameter to be passed on to BreakpointR, setting the maximum amount of background allowable for ww or cc calls.
#' @param type 'wc', 'ww', or both. Describes the input composite files (which are either Watson-Watson or Watson-Crick, or both).
#'
#' @return  A GRanges object of potential inversions for genoptyping with invertyper()
#'
#' @export
discover_possible_inversions <- function(composite_files=list(), windowsize=c(40,120,260), minReads=c(15,50,50), paired_reads=FALSE, numCPU=24, chromosomes=NULL, blacklist=NULL, background=0.2, type=c('wc','ww')){


stopifnot("The arguments minReads and windowsize should be parallel (have the same length)" = (length(minReads)==length(windowsize)))
stopifnot("This function should be used with two composite files named 'WC' and 'WW', for diploids (names(composite_files)). Also, a CC composite file won't work."  = (all(names(composite_files) %in% c('WC','WW')) & length(composite_files)==2) | all(type=='ww'))
stopifnot("This function should be used with a composite file named 'WW', for haploids (names(composite_files)). Also, a CC composite file won't work." = (all(names(composite_files) %in% c('WW')) & length(composite_files)==1) | any('wc' %in% type))

n <- names(composite_files)
message("begin granges conversion")
composite_files <- lapply(composite_files, galignment_to_granges, purpose='BreakpointR', paired_reads=paired_reads, region=NULL)
names(composite_files) <- n

		print(composite_files)

bpr <- list()
message("begin bpr loop")
for(i in c(1:length(windowsize))){


print(is.integer(minReads[i]))
print(is.integer(windowsize[i]))

message('replace binmethod size')
message(paste0("start bpr loop ",i))
bpr[[i]] <- breakpointr_for_invertyper(composite_files, plotspath=NULL,numCPU=numCPU, windowsize=as.integer(windowsize[i]), binMethod="size", minReads=as.integer(minReads[i]), background=background,maskRegions=blacklist,chromosomes=chromosomes)
message('remember to put back suppressMessages() above!!!')
gc()
}

message("end bpr loop")

if('wc' %in% type) {

possible_inversions <- sort(do.call('c',lapply(bpr, function(x)c(x$WC$counts[x$WC$counts$states!="wc"],x$WW$counts[x$WW$counts$states!="ww"]))))

} else {

possible_inversions <- sort(do.call('c',lapply(bpr, function(x)x$WW$counts[x$WW$counts$states!="ww"])))

}


message("time to return")

return(possible_inversions)

}


