
#' A stripped down parallel implementation of the core BreakpointR function, "runBreakpointr". 
#' Lots of this code is borrowed/altered from the BreakpointR package by David Porubsky, Ashley Sanders, and Aaron Taudt.
#' It also doesn't do such nice validation of arguments
#' Needs import foreach, doParallel, and maybe utils from read.table
#' Probably need to ad BreakpointR:: to all the BPR functions
#' We get a separate PDF for each file... maybe that could be improved?
#' qPDF could be used for that
#' plotspath should be a directory, not a filename
#' it uses single end GRanges with 2 frags paired.
#' For plotting, should also downsample beforehand to reduce processing time. Plus----sort the GRANGES! otherwise the chromosomes are all mixed up.




#' @param grangeslist
#' @param paired_reads
#' @param plotspath
#' @param numCPU
#' @param windowsize
#' @param binMethod
#' @param minReads
#' @param background
#' @param maskRegions
#' @param chromosomes
#'
#' @return  [...]
#'
#' @export
breakpointr_for_invertyper <- function(grangeslist=NULL, paired_reads=FALSE, plotspath=NULL,numCPU=4, windowsize=1000000, binMethod="size", minReads=50, background=0.2,maskRegions=NULL, chromosomes=NULL) {

	# so that R stops complaining about not being able to find dopar
	`%dopar%` <- foreach::`%dopar%`

	 runBreakpointr_for_InvertypeR <- function(granges, paired_reads=FALSE,plotspath=NULL,windowsize=1000000,binMethod="size",minReads=50,background=0.2,maskRegions=NULL, name='unknown', chromosomes=NULL){
	
		if(is.null(name)){
			name <- 'unknown'
		}

		tC <- tryCatch({
			breakpoints <- breakpointR::runBreakpointr(bamfile=granges, ID=name, pairedEndReads=paired_reads, pair2frgm=FALSE, chromosomes=chromosomes, windowsize=windowsize, binMethod=binMethod, background=background, 
						minReads=minReads,maskRegions=maskRegions)
			
			if(!is.null(plotspath)){
				if (!file.exists(plotspath)) {
    					dir.create(plotspath)
				}



				breakpoints$ID <- name
				breakpointR::plotBreakpoints(breakpoints, file=file.path(plotspath, paste0(name,".pdf")))
			}

			return(breakpoints)
		}, error = function(err) {
			stop(granges,'\n',err)
		})  
	
	}



	if (numCPU > 1) {

		## Parallelization ##
		message("Using ",numCPU," CPUs")
		cl <- parallel::makeCluster(numCPU)
		doParallel::registerDoParallel(cl)

		message("Finding breakpoints ...", appendLF=FALSE); ptm <- proc.time()
		temp <- foreach::foreach (index = c(1:length(grangeslist)), .packages='breakpointR', .final = function(x) setNames(x, names(grangeslist))) %dopar% {
			runBreakpointr_for_InvertypeR(granges = grangeslist[[index]], paired_reads=paired_reads, plotspath=plotspath,windowsize=windowsize,binMethod=binMethod,minReads=minReads,background=background,maskRegions=maskRegions, name=names(grangeslist)[index], chromosomes=chromosomes)
		}
    
		parallel::stopCluster(cl)
		time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

		} else {
			numCPU <- 1 # if to use only one CPU or CPU argument not defined
			temp <- list()

			for (index in c(1:length(grangeslist))) {
				temp[[index]] <- runBreakpointr_for_InvertypeR(granges = grangeslist[[index]], paired_reads=paired_reads,plotspath=plotspath,windowsize=windowsize, chromosomes=chromosomes,
						binMethod=binMethod,minReads=minReads,background=background,maskRegions=maskRegions, name=names(grangeslist)[index])
			}
			
			names(temp) <- names(grangeslist)
		}

	return(temp)
}


