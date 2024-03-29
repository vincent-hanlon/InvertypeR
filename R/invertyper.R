#' Genotypes inversions and adjusts their start and end coordinates.
#'
#' Requires WW and WC composite files (from create_composite_files() or in BAM format) and a GRanges object containing putative inversion coordinates
#'
#' We highly recommend using a hard_mask, and we provide one on the GitHub page for humans. For non-humans, regions that have anomalously high depth of coverage in many cells or that always have WC
#' (Watson-Crick) strand-state are suspect and should be added to a hard_mask by some custom analysis (for example, see the sup mat of https://doi.org/10.1186/s12864-021-07892-9).
#'
#' Priors should chosen to reflect the number of inversions you expect to find in an individual (up to ~100 if you use, for example, all publised inversion coordinates (2020)) and the number of unique
#  intervals in your BED file of putative inversions (many ways to calculate this, but 80% overlapping intervals could be considered the same, for example). If you are confident in the coordinates of
#' your putative inversions and you're just starting out, use adjust_method="raw". Otherwise, use adjust_method="deltas" or "all". The other options are specific to particular use cases.
#'
#' Note that we can phase inversions only because the WC composite file already has phased reads. This means that we know a 0|1 inversion on chr1 is on the same homolog as all 0|1 inversions
#' on chr1 in the sample, and that all chr1 1|0 inversions are on the other homolog. However, we don't know whether a 0|1 inversion on chr1 and a 0|1 inversion chr2 came from the same parent.
#' 0|1 inversions are distinguished from 1|0 inversions based on the strand switch in the WC composite file ( WC -> WW or WC -> CC).
#'
#' @param WW_reads A GenomicAlignments or GenomicAlignmentPairs object of reads from a Strand-seq WW composite file for an individual. Alternatively, the path to a BAM-formatted WW composite file.
#' @param WC_reads A GenomicAlignments or GenomicAlignmentPairs object of reads from a Strand-seq WC composite file for an individual. Alternatively, the path to a BAM-formatted WC composite file.
#' @param regions_to_genotype A GRanges object containing genomic intervals to be genotyped (putative inversions).
#' @param hard_mask A GRanges containing regions genomic intervals with poor-quality Strand-seq data. Reads that overlap these intervals will not be used. Highly recommended. Default "".
#' @param paired_reads Boolean: are the reads paired-end? Default TRUE.
#' @param confidence Posterior probability threshold above which you consider genotype calls to be reliable. Used to decide whether to keep adjusted inversions, as well as to identify low-confidence
#'   calls for adjust_method "low". Default 0.95.
#' @param prior Vector of three prior weights for inversion genotypes. For example, c("REF","HET","HOM") = c(0.9,0.05,0.05). Default c(0.33,0.33,0.33).
#' @param haploid_prior Vector of two prior weights for inversions on haploid chromosomes, like chrX and chrY in human males. For example, c("REF", "INV") = c(0.9,0.1). Default c(0.5,0.5).
#' @param chromosomes Vector of chromosome names to restrict the search for inversions.
#' @param haploid_chromosomes A vector of the names of chromosomes expected to be haploid (e.g., chrX and chrY in human males). Default NULL.
#' @param output_file Name of the file to write to. Default "inversions.txt".
#' @param adjust_method One of "raw", "merge", "deltas", "minimal", "low", or "all". Specifies which method to use to adjust the inversion coordinates 
#'   (start- and end-points). The adjustment routine ensures that adjusted inversions have the same genotype as they did before adjustment (if applicable), and that they 
#'   overlap at least one of the original unadjusted inversions in every cluster of overlapping events. If after adjustment (except with "raw") 
#'   overlapping inversions still remain, we merge the ones of the same genotype. If there are still overlapping inversions, we take the largest.
#'   "raw": No adjustment. "merge": Merge overlapping confident (i.e., posterior probability >= confidence for a genotype) inversions of the same genotype. "Deltas": adjust the endpoints of confident
#'   inversions based on the Strand-seq data. We calculate deltaW values for each read, which measure the change in read direction for nearby reads, and take the coordinates that correspond to the two
#'   highest values (i.e. the two spots with the greates change in read direction, which we hope will correspond to inversion breakpoints). This is done once for every set of overlapping inversions of the
#'   same genotype. "minimal": An experimental method. For confident overlapping inversions of the same genotype, take the interval common to all of them if one exists. "low": Much like "deltas", except
#'   that we only adjust inversions for which no confident genotype could be found (i.e., all prior probabilities < confidence). Sometimes this will "find" a confident inversion if the original coordinates
#'   were slightly off. "all": Does both "deltas" and "low". Default "all".
#' @param output_folder Directory for output files. Default "./".
#' @return A dataframe that associates each entry regions_to_genotype (putative inversions) with a genotype and posterior probability. Cols 1-3 are genomic coordinates, and cols 4-7 are read counts for the
#'   two composite files (WC and WW, or WC and CC if that's what you submitted) subset by Watson (W) and Crick (C) strands. Col 8 is the most probable genotype, and col 9 is the posterior probability
#'   associated with that genotype. Col 10 labels inversions with low read density. The start and end coordinates of such inversions should be checked by viewing Strand-seq data in a genome browser.
#'
#' @export
invertyper <- function(
    WW_reads,
    WC_reads,
    regions_to_genotype,
    hard_mask = NULL,
    paired_reads = TRUE,
    haploid_chromosomes = NULL,
    confidence = 0.95,
    prior = c(0.333, 0.333, 0.333),
    haploid_prior = c(0.5, 0.5),
    chromosomes = NULL,
    output_file = "inversions.txt",
    adjust_method = "all",
    output_folder = ".") {

    stopifnot(
        'Please choose a value for adjust_method. This controls how InvertypeR will attempt to improve the inversion coordinates you supply. Valid values are "raw","merge","deltas","minimal", "low", or "all".' =
            all(length(adjust_method) == 1 & adjust_method %in% c("raw", "merge", "deltas", "minimal", "low", "all"))
    )
    stopifnot(
        "Any haploid chromosomes you specify should be included as chromosomes too. This means including chrX and chrY as chromosomes AND as haploid_chromosomes for human males" =
            all(haploid_chromosomes %in% chromosomes) || is.null(chromosomes)
    )
    stopifnot("hard_mask should be a GRanges object, usually created from a BED file with import_bed()" = is.null(hard_mask) | class(hard_mask) == "GRanges")
    stopifnot("regions_to_genotype should be a GRanges object, usually created from a BED file with import_bed()" = is.null(regions_to_genotype) | class(regions_to_genotype) == "GRanges")

    haploid <- all(sort(haploid_chromosomes) == sort(chromosomes)) & !(is.null(chromosomes) & is.null(haploid_chromosomes))

    if((!all(sort(chromosomes) == sort(haploid_chromosomes)) | is.null(chromosomes) | is.null(haploid_chromosomes)) &
        all(prior == c(0.333, 0.333, 0.333))){
        warning(paste0("Using a default prior (c(",toString(prior),") for diploid chromosomes like the autosomes in humans) is not recommended. Consider what fraction of the putative inversions you wish to genotype are likely to have non-reference genotypes."))
    }

    if(!is.null(haploid_chromosomes) & all(haploid_prior == c(0.5, 0.5))){
        warning(paste0("Using a default haploid_prior (c(",toString(haploid_prior),") for haploid chromosomes like the sex chromosomes in human males) is not recommended. Consider what fraction of the putative inversions you wish to genotype are likely to have non-reference genotypes."))
    }

    # An irrelevant internal name change
    regions <- regions_to_genotype
    if (length(regions) == 0 || ((is(WW_reads, "GAlignmentPairs") || is(WW_reads, "GAlignments")) & length(WW_reads) < 5) || (((is(WC_reads, "GAlignmentPairs") || is(WC_reads, "GAlignments")) & length(WC_reads) < 5) & !haploid)) {
        warning("Non-empty composite files (min. 5 reads) and a non-empty GRanges for regions_to_genotype are required to genotype inversions")

        inversions <- data.frame(character(), integer(), integer(), integer(), integer(), integer(), integer(), character(), double())
        colnames(inversions) <- c("chr", "start", "end", "WW_counts_W", "WW_counts_C", "WC_counts_W", "WC_counts_C", "genotype", "probability")

        return(inversions)

    } else {
	
	if(!is.character(WW_reads)){
            found_chromosomes <- unique(sort(as.vector(S4Vectors::runValue(GenomicRanges::seqnames(WW_reads)))))
        } else {
            found_chromosomes <- unique(sort(as.vector(GenomicRanges::seqnames(regions_to_genotype))))
        }

        if (!is.null(chromosomes) & !is.character(WW_reads)) {
            if (!all(chromosomes %in% found_chromosomes)) {
                warning(paste0(
                    "Some chromosomes you provided (",
                    paste(chromosomes[!(chromosomes %in% found_chromosomes)], collapse = " "), ") don't have any reads in the composite files."
                ))
            }

            regions <- regions[as.vector(GenomicRanges::seqnames(regions)) %in% chromosomes]
        } else {
            chromosomes <- unique(sort(as.vector(GenomicRanges::seqnames(regions_to_genotype))))
        }

        if (is.character(WW_reads)) {
            file <- Rsamtools::BamFile(WW_reads)
            seqlengths_WW <- Rsamtools::scanBamHeader(file)$targets
            seqlengths_WC <- seqlengths_WW
        } else {
            seqlengths_WW <- GenomeInfoDb::seqlengths(WW_reads)
        }

        ptm <- startTimedMessage("\n       subsetting composite files and estimating background ...")

        # Takes reads from the BAM files that overlap the regions and stores them as GRanges
        WW_reads <- import_bam(WW_reads, chromosomes = chromosomes, paired_reads = paired_reads, hard_mask = hard_mask)
        if (!is.null(WC_reads)) {
            if (!is.character(WC_reads)) {
                seqlengths_WC <- GenomeInfoDb::seqlengths(WC_reads)
            }
            WC_reads <- import_bam(WC_reads, region = GenomicRanges::reduce(widen(
                granges = regions, seqlengths = seqlengths_WC,
                distance = 1e06
            ), min.gapwidth = 1000), paired_reads = paired_reads, hard_mask = hard_mask)
        } else {
            WC_reads <- WW_reads[0]
        }

        # Accurate background estimate, plus base strand state for the WW/CC file
        base <- WW_background(WW_reads, binsize = 1000000, paired_reads = paired_reads, chromosomes = found_chromosomes)

        WW_reads <- import_bam(WW_reads, region = GenomicRanges::reduce(widen(
            granges = regions, seqlengths = seqlengths_WW,
            distance = 1e06
        ), min.gapwidth = 1000), paired_reads = paired_reads, hard_mask = hard_mask)
        stopTimedMessage(ptm)
        # Genotyping the inversions
        msg <- paste0("       genotyping ", length(regions), " putative inversion(s) ...")

        ptm <- startTimedMessage(msg)
        inversions <- genotype_inversions(
            WW_reads = WW_reads, WC_reads = WC_reads, regions = regions, background = base[[1]], base_state = base[[2]], haploid_chromosomes = haploid_chromosomes,
            prior = prior, haploid_prior = haploid_prior
        )

        stopTimedMessage(ptm)

        # two methods under development for dealing with overlapping intervals
        if (adjust_method != "raw") {
            ptm <- startTimedMessage("       adjusting inversion coordinates ...")
            reads <- list(WW_reads, WC_reads)
            inversions <- adjust_intervals(
                inversions = inversions, reads = reads, confidence = confidence, base = base, haploid_chromosomes = haploid_chromosomes, prior = prior,
                haploid_prior = haploid_prior, adjust_method = adjust_method, paired_reads = paired_reads
            )
            inversions <- inversions[order(-inversions$probability), ]
            stopTimedMessage(ptm)
        }

        # Marking inversions with <0.001 reads per bp as being low-density, because they may span centromeres or chromosome ends/
        # The 0.001 figure is somewhat arbitrary, but it's adjusted based on the read count in the WW file of any submitted inversions.
	inversions[, 10] <- FALSE
        inversions[inversions[, 9] >= confidence & inversions[, 8] != 0 & inversions[, 8] != "0|0" & rowSums(inversions[, c(4:7)]) / (
            inversions[, 3] - inversions[, 2] + 1) < (0.001 * base[[3]] / 27734969), 10] <- TRUE

        colnames(inversions)[10] <- "low_read_density"

        if (!is.null(output_file)) {
            write.table(inversions, file.path(output_folder, output_file), quote = FALSE, sep = "\t", row.names = FALSE)
        }

        return(inversions)
    }
}
