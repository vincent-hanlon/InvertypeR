#' A stripped down parallel implementation of the core BreakpointR function, "runBreakpointr".
#'
#' Lots of this code is borrowed/altered from the BreakpointR package by David Porubsky, Ashley Sanders, and Aaron Taudt.
#' It doesn't do such nice validation of arguments etc., but it is sufficient for InvertypeR's purposes
#'
#'
#' @param grangeslist List of GRanges objects containing single-end reads (or the first read from paired reads)
#' @param plotspath File path where PDFs should be generated. Default NULL
#' @param numCPU See BreakpointR. Default 4.
#' @param windowsize See BreakpointR. Default 1000000.
#' @param binMethod See BreakpointR. Default "size".
#' @param minReads See BreakpointR. Default 50.
#' @param background See BreakpointR. Default 0.2.
#' @param maskRegions See BreakpointR. Default NULL.
#' @param chromosomes See BreakpointR. Default NULL.
#' @param parallelize_by_chromosome Set TRUE only if you wish to run different chromosomes on different threads, and then combine the results later. Use with care! 
#'     This will return a list() rather than a Breakpoint class object, it will ignore chromosomes with very few reads, and it will not combine plots or library metrics 
#'     nicely. It really is only compatible with discover_possible_inversions(). Default FALSE.
#' @return  A Breakpoints object from BreakpointR
#'
#' @export
breakpointr_for_invertyper <- function(grangeslist = NULL, plotspath = NULL, numCPU = 4, windowsize = 1000000, binMethod = "size", minReads = 50, background = 0.2, maskRegions = NULL, chromosomes = NULL, parallelize_by_chromosome = FALSE) {

    # so that R stops complaining about not being able to find dopar
    `%dopar%` <- foreach::`%dopar%`
    `%:%` <- foreach::`%:%`

    runBreakpointr_for_InvertypeR <- function(granges, paired_reads = FALSE, plotspath = NULL, windowsize = 1000000, binMethod = "size", minReads = 50, background = 0.2, maskRegions = NULL, name = "unknown", chromosomes = NULL) {
        if (is.null(name)) {
            name <- "unknown"
        }

        tC <- tryCatch(
            {
                breakpoints <- breakpointR::runBreakpointr(
                    bamfile = granges, ID = name, pairedEndReads = paired_reads, pair2frgm = FALSE, chromosomes = chromosomes, windowsize = windowsize, binMethod = binMethod, background = background,
                    minReads = minReads, maskRegions = maskRegions
                )

                if (!is.null(plotspath)) {
                    if (!file.exists(plotspath)) {
                        dir.create(plotspath)
                    }

                    breakpoints$ID <- name
                    breakpointR::plotBreakpoints(breakpoints, file = file.path(plotspath, paste0(name, ".pdf")))
                }

                return(breakpoints)
            },
            error = function(err) {
                stop(granges, "\n", err)
            }
        )
    }

    combine_Breakpoint_objects <- function(...){
        input <- list(...)
        seq <- unique(unlist(lapply(input, function(x)GenomeInfoDb::seqlevels(x$fragments))))
        fragments <- do.call('c', lapply(input, function(x){GenomeInfoDb::seqlevels(x$fragments) <- seq; x$fragments}))
        deltas <- do.call('c', lapply(input, function(x){GenomeInfoDb::seqlevels(x$deltas) <- seq; x$deltas}))
        breaks <- do.call('c', lapply(input, function(x){GenomeInfoDb::seqlevels(x$breaks) <- seq; x$breaks}))
        confint <- do.call('c', lapply(input, function(x){GenomeInfoDb::seqlevels(x$confint) <- seq; x$confint}))
        counts <- do.call('c', lapply(input, function(x){GenomeInfoDb::seqlevels(x$counts) <- seq; x$counts}))
        params <- input[[1]]$params
        lib.metrics <- input[[1]]$lib.metrics
        ID <- input[[1]]$ID
	
	bpr <- list(ID = ID, fragments = fragments, deltas = deltas, breaks = breaks, confint = confint, counts = counts, lib.metrics = lib.metrics, params = params)
        
        return(bpr)
    }

    if (numCPU > 1 & !parallelize_by_chromosome) {
        ## Parallelization ##
        message("Using ", numCPU, " CPUs")
        cl <- parallel::makeCluster(numCPU)
        doParallel::registerDoParallel(cl)

        message("Finding breakpoints ...", appendLF = FALSE)
        ptm <- proc.time()
        temp <- foreach::foreach(index = c(1:length(grangeslist)), .packages = "breakpointR", .final = function(x) setNames(x, names(grangeslist))) %dopar% {
            runBreakpointr_for_InvertypeR(granges = grangeslist[[index]], paired_reads = FALSE, plotspath = plotspath, windowsize = windowsize, binMethod = binMethod, minReads = minReads, background = background, maskRegions = maskRegions, name = names(grangeslist)[index], chromosomes = chromosomes)
        }

        parallel::stopCluster(cl)
        time <- proc.time() - ptm
        message(" ", round(time[3], 2), "s")
    } else if(numCPU > 1 & parallelize_by_chromosome) {
        # InvertypeR sometimes parallelizes by chromosome for extra speed
        if(is.null(chromosomes)){
            chromosomes <- sort(unlist(lapply(grangeslist, function(x)as.vector(GenomicRanges::seqnames(x)))))
            names(chromosomes) <- NULL
            chromosomes <- rle(chromosomes)
            indices <- which(chromosomes$lengths >= minReads)
            chromosomes <- chromosomes$values[indices]
        }
        cl <- parallel::makeCluster(numCPU)
        doParallel::registerDoParallel(cl)

        temp <- foreach::foreach(index = c(1:length(grangeslist)), .packages = "breakpointR") %:%
            foreach::foreach(chr=chromosomes, .packages = "breakpointR", .inorder=FALSE, .combine='combine_Breakpoint_objects', .multicombine=TRUE) %dopar% { 
                runBreakpointr_for_InvertypeR(granges = grangeslist[[index]][seqnames(grangeslist[[index]]) == chr], paired_reads = FALSE, plotspath = plotspath, windowsize = windowsize, binMethod = binMethod, minReads = minReads, background = background, maskRegions = maskRegions, name = names(grangeslist)[index], chromosomes = chr)
        }
        names(temp) <- names(grangeslist)
        parallel::stopCluster(cl)
    } else {
        numCPU <- 1 # if to use only one CPU or CPU argument not defined
        temp <- list()

        for (index in c(1:length(grangeslist))) {
            temp[[index]] <- runBreakpointr_for_InvertypeR(
                granges = grangeslist[[index]], paired_reads = FALSE, plotspath = plotspath, windowsize = windowsize, chromosomes = chromosomes,
                binMethod = binMethod, minReads = minReads, background = background, maskRegions = maskRegions, name = names(grangeslist)[index]
            )
        }

        names(temp) <- names(grangeslist)
    }

    return(temp)
}
