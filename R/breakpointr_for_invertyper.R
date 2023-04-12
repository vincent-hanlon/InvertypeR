#' A stripped down parallel implementation of the core BreakpointR function, "runBreakpointr".
#'
#' Lots of this code is borrowed/altered from the BreakpointR package by David Porubsky, Ashley Sanders, and Aaron Taudt.
#' It doesn't do such nice validation of arguments etc., but it is sufficient for InvertypeR's purposes
#'
#'
#' @param grangeslist List of GRanges objects containing single-end reads (or the first read from paired reads)
#' @param plotspath File path where PDFs should be generated. Default NULL
#' @param numCPU See BreakpointR
#' @param windowsize See BreakpointR
#' @param binMethod See BreakpointR
#' @param minReads See BreakpointR
#' @param background See BreakpointR
#' @param maskRegions See BreakpointR
#' @param chromosomes See BreakpointR
#'
#' @return  A Breakpoints object from BreakpointR
#'
#' @export
breakpointr_for_invertyper <- function(grangeslist = NULL, plotspath = NULL, numCPU = 4, windowsize = 1000000, binMethod = "size", minReads = 50, background = 0.2, maskRegions = NULL, chromosomes = NULL) {

    # so that R stops complaining about not being able to find dopar
    `%dopar%` <- foreach::`%dopar%`

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



    if (numCPU > 1) {
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
