#' Helps find reasonable priors for inversion genotyping
#'
#' With InvertypeR it is possible to genotype very long, messy lists of putative inversions from the literature, or from genomic features like inverted repeats
#' Sometimes it is difficult to choose a prior in these cases. What we do here is guess the number of inversions we think are present in the list of putative
#' inversions, and use that plus self-intersections of the inversion list to come up with a reasonable prior.
#'
#' This is not intended for cases where the inversion list does not have self-intersections----for example, if you are just genotyping a handful of well-characterized
#' inversions. Then you should create a prior based on known allele frequencies or something.
#'
#' This function could give ridiculous suggestions for priors if there are regions_to_genotype, for example. Use your judgement.
#'
#' @param regions_to_genotype A GRanges object of putative inversions. If you have a BED file, see import_bed()
#' @param expected_number_inversions An integer. If you had perfect data and genotyped an average individual using your regions_to_genotype, how many
#' unique inversions do you think there would be? For example, most humans have something like 100-200 paracentric inversions. If you think your inversion list
#' probably covers a lot of them, then maybe 100 is a good estimate. On the other hand, if your inversion list is a long shot and probably only a few of them are real,
#' you could consider a smaller number. Default 100.
#'
#' @return A list with the recommended prior and haploid_prior
#' @export
choose_priors <- function(regions_to_genotype, expected_number_inversions = 100) {

    hits <- GenomicRanges::findOverlaps(regions_to_genotype, regions_to_genotype)
    overlaps <- GenomicRanges::pintersect(regions_to_genotype[S4Vectors::queryHits(hits)], regions_to_genotype[S4Vectors::subjectHits(hits)])
    percentOverlap1 <- GenomicRanges::width(overlaps) / GenomicRanges::width(regions_to_genotype[S4Vectors::subjectHits(hits)])
    percentOverlap2 <- GenomicRanges::width(overlaps) / GenomicRanges::width(regions_to_genotype[S4Vectors::queryHits(hits)])
    intersections <- length(hits[percentOverlap1 > 0.8 & percentOverlap2 > 0.8])

    hom <- expected_number_inversions * intersections / (2 * length(regions_to_genotype)^2)
    het <- hom
    ref <- 1 - hom - het

    cat(paste0("Consider using the following priors:\n\tprior = c(", round(ref, 5), ", ", round(het, 5), ", ", round(hom, 5), ")\n\thaploid_prior = c(", round(ref, 5), ", ", round(het + hom, 5), ")\n"))

    return(invisible(list(c(ref, het, hom), c(ref, het + hom))))
}
