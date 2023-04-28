#' A pipeline-style implementation of all invertyper() functionality
#'
#' Requires Strand-seq BAM files, and preferably a VCF of heterozygous SNVs (can be called directly from the Strand-seq BAM files) and a BED file of regions to genotype
#'
#' In brief, this function:
#' 	(i) Creates phased Strand-seq composite files, which have higher coverage than any individual library.
#' 	(ii) Looks for strand-switches in the composite files that might be inversions, using BreakpointR
#' 	(iii) Genotypes the putative inversions that BreakpointR finds
#' 	(iv) Writes the inversions and relevant Strand-seq reads to UCSC Genome Browser files, if desired
#' 	(v) Genotypes any user-provided putative inversions. This typically helps find smaller inversions than those BreakpointR can detect
#' 	(vi) Writes those to UCSC Genome Browser files too
#' 	(vii) Returns dataframe(s) with inversion coordinates, genotypes, and posterior probabilities
#'
#' One tricky thing for human males is that it is necessary to set haploid_chromosomes=c('chrX','chrY')
#'
#' This was primarily set up for diploids, but haploids can be accomodated if haploid_chromosomes == chromosomes
#' @param regions_to_genotype A BED file containing genomic intervals to be genotyped (putative inversions). Default NULL.
#' @param hard_mask A BED file containing regions genomic intervals with poor-quality Strand-seq data. Reads that overlap these intervals will not be used. Highly recommended. Default NULL.
#' @param soft_mask A BED file containing regions with good Strand-seq data, but which interfere with composite file creation. These reads will appear in composite 
#'   files and inversion calls, but won't be used to identify regions with a given strand state. Typically these are large, obvious inversions or misorients, like the big chr8 
#'   inversion in humans. Initially, using the default NULL value is fine.
#' @param vcf A VCF file containing heterozygous SNVs for the sample. If there is no external source of SNVs, they can be called directly from the Strand-seq BAMs at a pinch
#' @param paired_reads Boolean: are the reads paired-end? Default TRUE.
#' @param confidence Posterior probability threshold above which you consider genotype calls to be reliable. Used to decide whether to keep adjusted inversions, as well as to identify low-confidence
#'   calls for adjust_method "low". Default 0.95.
#' @param prior Vector of three prior weights for inversion genotypes. For example, c("REF","HET","HOM") = c(0.9,0.05,0.05). Default c(0.33,0.33,0.33).
#' @param haploid_prior Vector of two prior weights for inversions on haploid chromosomes, like chrX and chrY in human males. For example, c("REF", "INV") = c(0.9,0.1). Default c(0.5,0.5).
#' @param chromosomes Vector of chromosome names to restrict the search for inversions.
#' @param haploid_chromosomes A vector of the names of chromosomes expected to be haploid (e.g., chrX and chrY in human males). Default NULL.
#' @param numCPU Integer. How many parallel threads to use, where possible?
#' @param input_folder Path to a directory containing Strand-seq BAM files
#' @param output_folder Path to a directory for output files
#' @param output_file Name of the inversion genotype file(s) to write to. Default "inversions.txt".
#' @param save_composite_files Should composite files be saved as RData objects?
#' @param write_browser_files Should inversions and associated reads be saved to UCSC Genome Browser files?
#' @param adjust_method One of "raw", "merge", "deltas", "minimal", "low", or "all". (default "all"). Specifies which method to use to adjust the inversion coordinates 
#'   (start- and end-points). The adjustment routine ensures that adjusted inversions have the same genotype as they 
#'   did before adjustment (if applicable), and that they overlap at least one of the original unadjusted inversions in every cluster of
#'   overlapping events. If after adjustment (except with "raw") overlapping inversions still remain, we merge the ones of the same genotype. If there are still overlapping inversions, we take the largest.
#'   "raw": No adjustment. "merge": Merge overlapping confident (i.e., posterior probability >= confidence for a genotype) inversions of the same genotype. "Deltas": adjust the endpoints of confident
#'   inversions based on the Strand-seq data. We calculate deltaW values for each read, which measure the change in read direction for nearby reads, and take the coordinates that correspond to the two
#'   highest values (i.e. the two spots with the greates change in read direction, which we hope will correspond to inversion breakpoints). This is done once for every set of overlapping inversions of the
#'   same genotype. "minimal": An experimental method. For confident overlapping inversions of the same genotype, take the interval common to all of them if one exists. "low": Much like "deltas", except
#'   that we only adjust inversions for which no confident genotype could be found (i.e., all prior probabilities < confidence). Sometimes this will "find" a confident inversion if the original coordinates
#'   were slightly off. "all": Does both "deltas" and "low".
#' @param discover_breakpointr_inversions Boolean. Should BreakpointR be used to find additional putative inversions not provided in regions_to_genotype? Default TRUE.
#' @param breakpointr_prior Vector of three prior weights for genotyping putative inversions found using BreakpointR. The default is ok for humans.
#' @param breakpointr_haploid_prior Vector of two prior weights. See breakpointr_prior
#' @param windowsize An integer vector, describing the binsizes for running BreakpointR, measured in # reads per bin. The default is ok for humans.
#' @param minReads An integer vector parallel to windowsize that gives the minimum number of reads allowed in a BreakpointR bin for deciding its strand state
#' @param background A BreakpointR parameter: how much background would you expect, at most, in a Watson-Watson region?
#'
#' @return  A dataframe or list of dataframe containing inversion coordinates, genotypes, posterior probabilities, and more.
#'
#' @export

invertyper_pipeline <- function(
    regions_to_genotype = NULL,
    prior = c(0.333, 0.333, 0.333),
    haploid_prior = c(0.5, 0.5),
    adjust_method = "all",
    input_folder = "./",
    output_folder = "./",
    haploid_chromosomes = NULL,
    vcf = NULL,
    paired_reads = TRUE,
    confidence = 0.95,
    hard_mask = NULL,
    soft_mask = NULL,
    chromosomes = NULL,
    numCPU = 4,
    save_composite_files = FALSE,
    write_browser_files = FALSE,
    discover_breakpointr_inversions = TRUE,
    breakpointr_prior = c(0.9, 0.05, 0.05),
    breakpointr_haploid_prior = c(0.9, 0.1),
    windowsize = c(40, 120, 360),
    minReads = c(15, 50, 50),
    background = 0.2,
    output_file = "inversions.txt") {
 
    stopifnot(
        "There is nothing to genotype: either provide a list of putative inversions (regions_to_genotype) or set discover_breakpointr_inversions=TRUE" =
            !is.null(regions_to_genotype) | discover_breakpointr_inversions
    )
    stopifnot(
        "The confidence threshold for posterior probabilties must be between 0 and 1 (and it should be of type numeric). In general, the arbitrary convention is to choose 0.95." =
            (as.numeric(confidence) < 1 & as.numeric(confidence) > 0 & is.numeric(confidence))
    )
    stopifnot(
        'Please choose a value for adjust_method. This controls how InvertypeR will attempt to improve the inversion coordinates you supply. Valid values are "raw","merge","deltas","minimal", "low", or "all".' =
            all(length(adjust_method) == 1 & adjust_method %in% c("raw", "merge", "deltas", "minimal", "low", "all")) | is.null(regions_to_genotype)
    )
    stopifnot("The arguments minReads and windowsize should be parallel (have the same length)" = (length(minReads) == length(windowsize)))
    stopifnot("The regions_to_genotype file cannot be found" = is.null(regions_to_genotype) || file.exists(regions_to_genotype))
    stopifnot("The hard_mask file cannot be found" = is.null(hard_mask) || file.exists(hard_mask))
    stopifnot("The soft_mask file cannot be found" = is.null(soft_mask) || file.exists(soft_mask))
    stopifnot("The vcf file cannot be found" = is.null(vcf) || file.exists(vcf))
    stopifnot("The input_folder cannot be found" = file.exists(input_folder))
    stopifnot(
        "Any haploid chromosomes you specify should be include as chromosomes too. This means including chrX and chrY as chromosomes AND as haploid_chromosomes for human males" =
            all(haploid_chromosomes %in% chromosomes) || is.null(chromosomes)
    )

    if (!is.null(regions_to_genotype) & (all(prior == c(0.333, 0.333, 0.333)) & (!all(sort(chromosomes) == sort(haploid_chromosomes)) | is.null(haploid_chromosomes)) |
        (!is.null(haploid_chromosomes) & all(haploid_prior == c(0.5, 0.5))))) {
        warning("Using the default priors (prior for homogametic diploids (e.g., human females), haploid_prior for haploids, or both for heterogametic dipoids (e.g., human males)) is not recommended. Consider what fraction of the putative inversions you wish to genotype are likely to have non-reference genotypes.")
    }

    if (!is.null(regions_to_genotype)) {
        regions_to_genotype <- import_bed(regions_to_genotype)
    }

    if (!is.null(hard_mask)) {
        hard_mask <- import_bed(hard_mask)
    }

    if (!is.null(soft_mask)) {
        soft_mask <- import_bed(soft_mask)
    }

    type <- c("wc", "ww")

    if (!is.null(haploid_chromosomes) & !is.null(chromosomes) & all(sort(chromosomes) == sort(haploid_chromosomes))) {
        type <- "ww"
    }

    composite_files <- create_composite_files(
        input_folder = input_folder, type = type, numCPU = numCPU, vcf = vcf, paired_reads = paired_reads, hard_mask = hard_mask,
        soft_mask=soft_mask, output_folder = output_folder, save_composite_files = save_composite_files, chromosomes = chromosomes
    )

    for (i in composite_files) {
        i <- drop_seq_qual(i, paired_reads = paired_reads)
    }

    if (discover_breakpointr_inversions) {
        possible_inversions <- discover_possible_inversions(
            composite_files = composite_files, windowsize = windowsize, minReads = minReads, paired_reads = paired_reads,
            numCPU = numCPU, chromosomes = chromosomes, hard_mask = hard_mask, background = background, type = type
        )

        if (length(composite_files) == 1) {
            composite_files[2] <- list(NULL)
            names(composite_files)[2] <- "WC"
        }

        breakpointr_inversions <- invertyper(
            WW_reads = composite_files$WW, WC_reads = composite_files$WC, regions_to_genotype = possible_inversions, hard_mask = hard_mask, paired_reads = paired_reads, 
            haploid_chromosomes = haploid_chromosomes, confidence = confidence, prior = breakpointr_prior,
            haploid_prior = breakpointr_haploid_prior, output_file = paste0("breakpointr_", output_file), adjust_method = c("merge"), output_folder = output_folder
        )

        if (write_browser_files) {
            write_UCSC_browser_files(
                inversions = breakpointr_inversions, WW_reads = composite_files$WW, WC_reads = composite_files$WC,
                confidence = confidence, paired_reads = paired_reads, output_folder = output_folder, prefix = "breakpointr_", type = type
            )
        }
    }

    if (!is.null(regions_to_genotype)) {
        inversions <- invertyper(
            WW_reads = composite_files$WW, WC_reads = composite_files$WC, regions_to_genotype = regions_to_genotype, hard_mask = hard_mask, paired_reads = paired_reads, 
            haploid_chromosomes = haploid_chromosomes, confidence = confidence, prior = prior, haploid_prior = haploid_prior,
            output_file = output_file, adjust_method = adjust_method, output_folder = output_folder
        )

        if (write_browser_files) {
            write_UCSC_browser_files(
                inversions = inversions, WW_reads = composite_files$WW, WC_reads = composite_files$WC, confidence = confidence,
                paired_reads = paired_reads, output_folder = output_folder, prefix = "", type = type
            )
        }
    }

    if (!is.null(regions_to_genotype) & discover_breakpointr_inversions) {
        return(list(inversions, breakpointr_inversions))
    } else if (!is.null(regions_to_genotype)) {
        return(inversions)
    } else {
        return(breakpointr_inversions)
    }
}
