#' Discovering putative inversions as strand-switches with BreakpointR
#'
#' @param composite_files A list of 1 or 2 composite files, named WC and WW (or just WW for haploids). Default list().
#' @param windowsize A vector of integers. Each integer will be used the set the number of reads in a bin for a BreakpointR run. Default c(40, 120, 260).
#' @param minReads A vector of integers parallel to windowsize, giving the minimum number of reads in a bin for it to be considered/genotyped. Default c(15, 50, 50).
#' @param paired_reads Boolean. Paired-end reads? Default FALSE.
#' @param numCPU Integer. How many threads to use? Default 4.
#' @param chromosomes A character vector of chromosome names to examine
#' @param hard_mask A GRanges object containing intervals with suspected poor-quality Strand-seq data. 
#' @param background A parameter to be passed on to BreakpointR, setting the maximum amount of background allowable for ww or cc calls. Default 0.2.
#' @param type 'wc', 'ww', or both. Describes the input composite files (which are either Watson-Watson or Watson-Crick, or both). Default c("wc", "ww").
#'
#' @return  A GRanges object of potential inversions for genoptyping with invertyper()
#'
#' @export
discover_possible_inversions <- function(
    composite_files = list(),
    windowsize = c(40, 120, 260),
    minReads = c(15, 50, 50),
    paired_reads = FALSE,
    numCPU = 4,
    chromosomes = NULL,
    hard_mask = NULL,
    background = 0.2,
    type = c("wc", "ww")) {

    stopifnot("The arguments minReads and windowsize should be parallel (have the same length)" = (length(minReads) == length(windowsize)))
    stopifnot("hard_mask should be a GRanges object, usually created from a BED file with import_bed()" = is.null(hard_mask) | class(hard_mask) == "GRanges")
    stopifnot(
        "This function should be used with two composite files named 'WC' and 'WW', for diploids (names(composite_files)). Also, a CC composite file won't work." =
            (all(names(composite_files) %in% c("WC", "WW")) & length(composite_files) == 2) | all(type == "ww")
    )
    stopifnot(
        "This function should be used with a composite file named 'WW', for haploids (names(composite_files)). Also, a CC composite file won't work." =
            (all(names(composite_files) %in% c("WW")) & length(composite_files) == 1) | any("wc" %in% type)
    )

    n <- names(composite_files)
    composite_files <- lapply(composite_files, galignment_to_granges, purpose = "BreakpointR", paired_reads = paired_reads, region = NULL)
    names(composite_files) <- n

    bpr <- list()
    start <- ""

    for (i in c(1:length(windowsize))) {
        if (i > 3) {
            suffix <- "th"
        } else if (i == 3) {
            suffix <- "rd"
        } else if (i == 2) {
            suffix <- "nd"
        } else if (i == 1) {
            suffix <- "st"
            start <- "\n"
        }

        msg <- paste0(start, "       running BreakpointR to discover inversions for the ", i, suffix, " time ...")
        ptm <- startTimedMessage(msg)

        bpr[[i]] <- suppressMessages(breakpointr_for_invertyper(composite_files,
            plotspath = NULL, numCPU = numCPU, windowsize = as.integer(windowsize[i]), binMethod = "reads",
            minReads = as.integer(minReads[i]), background = background, maskRegions = hard_mask, chromosomes = chromosomes,
            parallelize_by_chromosome = TRUE
        ))
        invisible(gc())
        stopTimedMessage(ptm)

        start <- ""
    }

    if ("wc" %in% type) {
        possible_inversions <- sort(do.call("c", lapply(bpr, function(x) c(x$WC$counts[x$WC$counts$states != "wc"], x$WW$counts[x$WW$counts$states != "ww"]))))
    } else {
        possible_inversions <- sort(do.call("c", lapply(bpr, function(x) x$WW$counts[x$WW$counts$states != "ww"])))
    }

    return(possible_inversions)
}
