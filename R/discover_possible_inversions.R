#' @param composite_files
#' @param windowsize
#' @param minReads
#' @param paired_reads
#' @param numCPU
#' @param chromosomes
#' @param blacklist
#' @param background
#'
#' @return  [...]
#'
#' @export
discover_possible_inversions <- function(composite_files=list(), windowsize=c(40,120,260), minReads=c(15,50,50), paired_reads=FALSE, numCPU=24, chromosomes=NULL, blacklist=NULL, background=0.2){


stopifnot("The arguments minReads and windowsize should be parallel (have the same length)" = (length(minReads)==length(windowsize)))
stopifnot("This function can only be used with a pair of composite files, named 'WC' and 'WW' (names(composite_files)). Also, a CC composite file won't work."  = (all(names(composite_files) %in% c('WC','WW')) & length(composite_files)==2))

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


message(paste0("start bpr loop ",i))
bpr[[i]] <- breakpointr_for_invertyper(composite_files, paired_reads=FALSE, plotspath=NULL,numCPU=numCPU, windowsize=as.integer(windowsize[i]), binMethod="reads", minReads=as.integer(minReads[i]), background=background,maskRegions=blacklist,chromosomes=chromosomes)
message('remember to put back suppressMessages() above!!!')
gc()
}

message("end bpr loop")
possible_inversions <- sort(do.call('c',lapply(bpr, function(x)c(x$WC$counts[x$WC$counts$states!="wc"],x$WW$counts[x$WW$counts$states!="ww"]))))

message("time to return")
return(possible_inversions)

}


