#' Creates Strand-seq composite files for inversion genotyping
#'
#' This replaces BASH scripts in a previous InvertypeR version.
#'
#' In brief, this function merges Strand-seq libraries into two composite files. The Watson-Watson (WW, or CC) composite file contains reads from Watson-Watson regions,
#' plus reads from Crick-Crick regions that have been reversed before merging. The Watson-Crick composite file combines reads from Watson-Crick regions and 'Crick-Watson'
#' regions. Basically, this last bit is a haplotype aware merging step that accounts for the fact that half of the time, a Watson-Crick/Crick-Watson region has reads from
#' homolog 1 align in the reverse orientation and reads from homolog 2 align in the forward orientation, while half of the time the orientations are switched.
#'
#' The composite files can be saved to an RData object for later use. However, since the composite files are GenomicAlignments or GenomicAlignmentPairs objects, it should
#' also be possible to save them as samtools-compatible BAM files using rtracklayer::export().
#'
#' @param input_folder The path to a directory containing good-quality Strand-seq BAM files
#' @param type Either 'ww', or 'wc', or both. The 'wc' file is only available for diploids for which heterozygous SNPs are available, but it greatly improves inversion calls
#' @param numCPU An integer; the number of threads to use
#' @param vcf The path to a VCF file containing heterozygous SNPs for an individual. External SNPs (e.g. from WGS data) are best, but they can also usually be called
#' directly from the Strand-seq data using BBTOOLS callvariants.sh (see README on GitHub).
#' @param paired_reads Boolean. Are the reads paired-end?
#' @param hard_mask A GRanges object containing regions that are thought to contain unreliable Strand-seq data. Highly recommended.
#' @param soft_mask A GRanges object containing regions with good Strand-seq data, but which interfere with composite file creation. These reads will appear in composite files and inversion calls, but won't be used to identify regions with a given strand state. Typically these are large, obvious inversions or misorients, like the big chr8 inversion in humans. Initially, using the default NULL value is fine.
#' @param output_folder The name of a folder to which output files should be written.
#' @param save_composite_files Boolean. Should composite files be saved as RData objects, or not at all?
#' @param chromosomes A character vector of chromosome names for which composite files should be created.
#'
#' @return  A list containing 1 or 2 composite files in GenomicAlignments format.
#'
#' @export
create_composite_files <- function(
    input_folder = "./",
    type = c("wc", "ww"),
    numCPU = 4,
    vcf = NULL,
    paired_reads = TRUE,
    hard_mask = NULL,
    soft_mask = NULL,
    output_folder = "./",
    save_composite_files = FALSE,
    chromosomes = NULL) {

    stopifnot(
        "The type argument should be 'wc' for a Watson-Crick composite file, 'ww' for a Watson-Watson composite file, or c('ww','wc') for both." =
            (!any(!type %in% c("wc", "ww")))
    )
    stopifnot("The path to a VCF file is required for a Watson-Crick ('wc') composite file." = !(is.null(vcf) & "wc" %in% type))

    if (is.null(hard_mask)) {
        warning("A hard_mask is highly recommended both for creating composite files and for genotyping inversions.")
    }

    if (!file.exists(output_folder)) {
        dir.create(output_folder)
    }

    ptm <- startTimedMessage("\n       loading BAM files ...")
    bamlist <- list.files(input_folder, pattern = "\\.bam$")
    stopifnot("No BAM files found in input_folder for composite file creation." = length(bamlist) > 0)

    cl <- suppressMessages(parallel::makeCluster(numCPU))
    galignmentslist <- parallel::parLapply(cl, bamlist, import_bam, hard_mask = hard_mask, chromosomes = chromosomes, paired_reads = paired_reads)

    names(galignmentslist) <- bamlist
    found_chromosomes <- unique(sort(do.call("c", lapply(galignmentslist, function(x) unique(sort(as.vector(S4Vectors::runValue(GenomicRanges::seqnames(x)))))))))

    if (!all(chromosomes %in% found_chromosomes)) {
        warning(paste0("Some chromosomes you provided (", paste(chromosomes[!(chromosomes %in% found_chromosomes)], collapse = " "), ") don't have any reads in the BAM files."))
    }

    grangeslist <- parallel::parLapply(cl, galignmentslist, galignment_to_granges, purpose = "BreakpointR", paired_reads = paired_reads, pair2frgm = FALSE)
    names(grangeslist) <- names(galignmentslist)
    parallel::stopCluster(cl)

    stopTimedMessage(ptm)
    ptm <- startTimedMessage("       recording the strand states of genomic regions ...")

    mask <- NULL
    if (!is.null(soft_mask) & !is.null(hard_mask)) {
        mask <- GenomicRanges::reduce(GenomicRanges::sort(c(hard_mask, soft_mask)))
    } else if (!is.null(soft_mask)) {
        mask <- soft_mask
    } else if (!is.null(hard_mask)) {
        mask <- hard_mask
    }

    bpr <- suppressMessages(breakpointr_for_invertyper(grangeslist,
        numCPU = numCPU, windowsize = 20000000, binMethod = "size", minReads = 50, background = 0.2, maskRegions = mask,
        chromosomes = chromosomes
    ))

    rm("grangeslist")
    invisible(gc())

    WWregions <- sort(find_regions_with_strand_state(bpr, states = c("ww", "cc"), region_size = 100000))
    WWregions$Cs <- NULL
    WWregions$Ws <- NULL
    WWregions$state <- NULL
    GenomeInfoDb::seqlevels(WWregions) <- GenomeInfoDb::seqlevels(galignmentslist[[1]])

    if ("wc" %in% type) {
        WCregions <- sort(find_regions_with_strand_state(bpr, states = c("wc"), region_size = 100000))
        WCregions$Cs <- NULL
        WCregions$Ws <- NULL
        WCregions$state <- NULL
        GenomeInfoDb::seqlevels(WCregions) <- GenomeInfoDb::seqlevels(galignmentslist[[1]])
        rm("bpr")
        invisible(gc())

        all_phased_WCregions <- strandPhaseR_for_invertyper(
            numCPU = numCPU, positions = vcf, WCregions = WCregions, chromosomes = chromosomes,
            paired_reads = paired_reads, num.iterations = 3, galignmentslist = galignmentslist
        )
    } else {
        rm("bpr")
        invisible(gc())
    }
    stopTimedMessage(ptm)
    ptm <- startTimedMessage("       combining reads into composite files ...")

    if (.Platform$OS.type == "windows") {
        mc.cores <- 1
    } else {
        mc.cores <- numCPU
    }

    # If this is parallelized with parLapply instead it might work on Windows, but unfortunately parLapply can't subset GenomicRanges like granges[vector] for some weird reason.
    CC <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state,
        galignmentslist = galignmentslist, regions = WWregions,
        paired_reads = paired_reads, states = "cc", flip_reads = TRUE, mc.cores = mc.cores
    )
    WW <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state,
        galignmentslist = galignmentslist, regions = WWregions,
        paired_reads = paired_reads, states = "ww", flip_reads = FALSE, mc.cores = mc.cores
    )

    if ("wc" %in% type) {
        CW <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state,
            galignmentslist = galignmentslist,
            regions = all_phased_WCregions, paired_reads = paired_reads, states = "cw", flip_reads = TRUE, mc.cores = mc.cores
        )
        WC <- parallel::mclapply(names(galignmentslist), extract_reads_by_region_and_state,
            galignmentslist = galignmentslist,
            regions = all_phased_WCregions, paired_reads = paired_reads, states = "wc", flip_reads = FALSE, mc.cores = mc.cores
        )
    }

    rm("galignmentslist")
    invisible(gc())

    WW_composite_file <- do.call("c", c(CC, WW))

    to_plot <- list(galignment_to_granges(WW_composite_file[sort(sample(length(WW_composite_file), min(length(WW_composite_file), 500000)))],
        purpose = "BreakpointR", paired_reads = paired_reads, pair2frgm = FALSE
    ))
    names(to_plot)[1] <- "WW_composite_file"

    composite <- list(WW_composite_file)
    names(composite) <- c("WW")

    if ("wc" %in% type) {
        WC_composite_file <- do.call("c", c(WC, CW))
        to_plot <- append(to_plot, galignment_to_granges(WC_composite_file[sort(sample(length(WC_composite_file), min(length(WC_composite_file), 500000)))],
            purpose = "BreakpointR", pair2frgm = FALSE, paired_reads = paired_reads
        ))

        names(to_plot)[2] <- "WC_composite_file"

        composite <- append(composite, WC_composite_file)

        names(composite) <- c("WW", "WC")
    }

    invisible(suppressMessages(breakpointr_for_invertyper(to_plot,
        plotspath = output_folder, numCPU = numCPU, windowsize = 1000000,
        binMethod = "size", minReads = 50, background = 0.2, maskRegions = NULL, chromosomes = chromosomes
    )))
    stopTimedMessage(ptm)

    if (save_composite_files) {
        ptm <- startTimedMessage("       saving composite files ...")
        if ("ww" %in% type) {
            save(WW_composite_file, file = file.path(output_folder, "WW_composite_file.RData"))
        }

        if ("wc" %in% type) {
            save(WC_composite_file, file = file.path(output_folder, "WC_composite_file.RData"))
        }
        stopTimedMessage(ptm)
    }

    return(composite)
}
