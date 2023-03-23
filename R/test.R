#' test function
#'
#' @export
test <- function(granges) {


message(paste(c(deparse(sys.calls()[[sys.nframe()-1]]),deparse(substitute(granges)) ,as.character(length(GenomicRanges::trim(granges))==length(granges)))))
}


