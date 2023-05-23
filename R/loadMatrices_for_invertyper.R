#' Copied from StrandPhaseR so that it accesses InvertypeR's version of bamregion2granges
#' Vincent Hanlon changed @inheritParams, the bamregions function, some added package:: tags, and some improvements to memory efficiency
#'
#' This funcion will read in partial single cell haplotypes into
#' two parallel matrices separately for Watson and Crick reads
#'
#' @inheritParams phaseChromosome_for_invertyper
#' @importFrom GenomicAlignments pileLettersAt
#' @importFrom Biostrings alphabetFrequency DNAStringSet
#'
#' @author David Porubsky
loadMatrices_for_invertyper <- function(input_folder = NULL, positions = NULL, WCregions = NULL, pairedEndReads = FALSE, min.mapq = 10, min.baseq = 30, galignmentslist=galignmentslist) {

    ## function to collapse list of letters
    collapse.str <- function(x) {
        if (length(x) > 1) {
            paste(x, collapse = "")
        } else if (length(x) == 1) {
            as.character(x)
        } else {
            x <- ""
        }
    }

    message(" Loading data for ", length(WCregions), " WCregions", appendLF = F)
    ptm <- proc.time()

    file.list <- unique(WCregions$filename)

    # VH: Making loops into functions to save a bit of memory
    covered_snvs_file <- function(file, input_folder, positions, pairedEndReads, min.mapq) {
        # message("Processing ", file)

        bam.reads <- bamregion2GRanges_for_invertyper(file, region = positions, pairedEndReads = pairedEndReads, min.mapq = min.mapq, galignmentslist=galignmentslist)
        seqlengths(WCregions) <- seqlengths(bam.reads)

        regions <- WCregions[WCregions$filename == file]

        ## process each WC region at a time

        covered <- lapply(c(1:length(regions)), covered_snvs_region, regions = regions, input_folder = input_folder, positions = positions, bam.reads = bam.reads)
        covered <- do.call(Map, c(c, covered))
        return(covered)
    }

    # VH: Making loops into functions again...
    covered_snvs_region <- function(i, regions, input_folder, positions, bam.reads) {
        genomic.pos <- vector()
        watson.bases <- list()
        crick.bases <- list()
        watson.quals <- list()
        crick.quals <- list()

        region <- regions[i]
        # message("Working on ", region)

        ## make unique ID per WC region
        filename <- region$filename
        region.ID <- as.character(region)
        filename.ID <- paste(filename, region.ID, sep = "__")

        bamfile <- file.path(input_folder, filename)

        ## get snv position overlaping with given WC region
        mask <- GenomicRanges::findOverlaps(region, positions)
        region.snvs <- positions[S4Vectors::subjectHits(mask)]
        region.snvs <- sort(region.snvs)

        ## make sure that there are no duplicated snv positions in region.snv list
        region.snvs <- region.snvs[!duplicated(GenomicRanges::start(region.snvs))]

        ## Get data from the list of GRanges
        bam.overlaps <- GenomicRanges::findOverlaps(region, bam.reads)
        data <- bam.reads[S4Vectors::subjectHits(bam.overlaps)]
        seqlevels(data) <- seqlevels(region)

        ## Split reads by directionality
        crick <- data[strand(data) == "+"]
        watson <- data[strand(data) == "-"]

        ## extract read sequences for watson and crick reads
        crick.seq <- GenomicRanges::mcols(crick)$seq
        watson.seq <- GenomicRanges::mcols(watson)$seq

        ## extract base qualities for watson and crick reads
        crick.qual <- GenomicRanges::mcols(crick)$qual
        watson.qual <- GenomicRanges::mcols(watson)$qual

        ## get piles of bases at each variable position
        piles.crick <- GenomicAlignments::pileLettersAt(crick.seq, GenomicRanges::seqnames(crick), GenomicRanges::start(crick), GenomicRanges::mcols(crick)$cigar, region.snvs)
        piles.watson <- GenomicAlignments::pileLettersAt(watson.seq, GenomicRanges::seqnames(watson), GenomicRanges::start(watson), GenomicRanges::mcols(watson)$cigar, region.snvs)

        ## get piles of qualities at each variable position
        quals.crick <- GenomicAlignments::pileLettersAt(crick.qual, GenomicRanges::seqnames(crick), GenomicRanges::start(crick), GenomicRanges::mcols(crick)$cigar, region.snvs)
        quals.watson <- GenomicAlignments::pileLettersAt(watson.qual, GenomicRanges::seqnames(watson), GenomicRanges::start(watson), GenomicRanges::mcols(watson)$cigar, region.snvs)

        ## Filter bases based on base quality
        df.pilesCrick <- as(piles.crick, "data.frame")
        df.pilesWatson <- as(piles.watson, "data.frame")
        df.qualsCrick <- as(quals.crick, "data.frame")
        df.qualsWatson <- as(quals.watson, "data.frame")

        # bases.crick <- strsplit(df.pilesCrick$x, "")
        # bases.watson <- strsplit(df.pilesWatson$x, "")
        bases.crick <- strsplit(df.pilesCrick[, 1], "")
        bases.watson <- strsplit(df.pilesWatson[, 1], "")

        ## translate raw base qualities into a number (for Sanger 33)
        # quals.crick <- sapply(df.qualsCrick$x, function(x) as.numeric(charToRaw(x))-33)
        # quals.watson <- sapply(df.qualsWatson$x, function(x) as.numeric(charToRaw(x))-33)
        quals.crick <- sapply(df.qualsCrick[, 1], function(x) as.numeric(charToRaw(x)) - 33)
        quals.watson <- sapply(df.qualsWatson[, 1], function(x) as.numeric(charToRaw(x)) - 33)

        ## TODO: change mean base quality values to calculate join probability for overlapping reads
        filtbases.crick <- mapply(function(X, Y) {
            X[Y >= min.baseq]
        }, X = bases.crick, Y = quals.crick)
        filtquals.crick <- sapply(quals.crick, function(x) round(mean(x[x >= min.baseq])))
        filtbases.watson <- mapply(function(X, Y) {
            X[Y >= min.baseq]
        }, X = bases.watson, Y = quals.watson)
        filtquals.watson <- sapply(quals.watson, function(x) round(mean(x[x >= min.baseq])))

        filtbases.crick <- lapply(filtbases.crick, collapse.str)
        filtbases.watson <- lapply(filtbases.watson, collapse.str)
        filtbases.crick <- Biostrings::DNAStringSet(unlist(filtbases.crick))
        filtbases.watson <- Biostrings::DNAStringSet(unlist(filtbases.watson))

        ## calculate base frequency at each variable position
        crickBaseFreq.m <- Biostrings::alphabetFrequency(filtbases.crick, baseOnly = T) # get base frequency matrix
        watsonBaseFreq.m <- Biostrings::alphabetFrequency(filtbases.watson, baseOnly = T)

        ## delete snv position with more than one alternative allele
        crick.mask <- rowSums(crickBaseFreq.m > 0)
        watson.mask <- rowSums(watsonBaseFreq.m > 0)
        crickBaseFreq.m[crick.mask > 1, ] <- 0 # positions with more than one alternative allele set to zero
        watsonBaseFreq.m[watson.mask > 1, ] <- 0

        ## filter region with only one snv coveraged either on W or C reads
        crickBaseCov <- rowSums(crickBaseFreq.m)
        watsonBaseCov <- rowSums(watsonBaseFreq.m)
        if (length(crickBaseCov[crickBaseCov > 0]) < 2 | length(watsonBaseCov[watsonBaseCov > 0]) < 2) {
            return(list(list(),list(),list(),list(),list()))
        }

        ## get covered positions
        crickBaseFreq.collapsed <- which(crickBaseFreq.m > 0, arr.ind = T) # find postions of covered SNVs and column index for A,C,G,T as 1,2,3,4
        crickBaseFreq.collapsed.srt <- crickBaseFreq.collapsed[order(crickBaseFreq.collapsed[, 1]), ] # sort by SNV postion
        crickBaseFreq.collapsed.srt <- crickBaseFreq.collapsed.srt[crickBaseFreq.collapsed.srt[, 2] <= 4, ] # remove unusual bases having values other than 1,2,3,4 (N's)
        crick.pos <- crickBaseFreq.collapsed.srt[, 1]

        watsonBaseFreq.collapsed <- which(watsonBaseFreq.m > 0, arr.ind = T)
        watsonBaseFreq.collapsed.srt <- watsonBaseFreq.collapsed[order(watsonBaseFreq.collapsed[, 1]), ]
        watsonBaseFreq.collapsed.srt <- watsonBaseFreq.collapsed.srt[watsonBaseFreq.collapsed.srt[, 2] <= 4, ]
        watson.pos <- watsonBaseFreq.collapsed.srt[, 1]

        ## transform obtained matrix into a vector with covered snv positions as values and base IDs as names
        crickBases.v <- GenomicRanges::start(region.snvs[crick.pos])
        crickBases.df <- as.data.frame(crickBaseFreq.collapsed.srt)
        names(crickBases.v) <- crickBases.df[crickBases.df$row %in% crick.pos, ]$col

        # watsonBases.v <- watson.pos
        watsonBases.v <- GenomicRanges::start(region.snvs[watson.pos])
        watsonBases.df <- as.data.frame(watsonBaseFreq.collapsed.srt)
        names(watsonBases.v) <- watsonBases.df[watsonBases.df$row %in% watson.pos, ]$col

        ## transform base qualities into a vector
        crickQuals.v <- GenomicRanges::start(region.snvs[crick.pos])
        watsonQuals.v <- GenomicRanges::start(region.snvs[watson.pos])
        names(crickQuals.v) <- filtquals.crick[crick.pos]
        names(watsonQuals.v) <- filtquals.watson[watson.pos]

        ## store genomic positions of covered SNVs
        genomic.pos <- c(genomic.pos, crickBases.v, watsonBases.v)

        crick.bases[[filename.ID]] <- crickBases.v
        watson.bases[[filename.ID]] <- watsonBases.v
        crick.quals[[filename.ID]] <- crickQuals.v
        watson.quals[[filename.ID]] <- watsonQuals.v

        return(list(genomic.pos, crick.bases, watson.bases, crick.quals, watson.quals))
    }



    covered <- lapply(file.list, covered_snvs_file, input_folder = input_folder, positions = positions, pairedEndReads = pairedEndReads, min.mapq = min.mapq)


    covered <- do.call(Map, c(c, covered))

    genomic.pos <- unique(covered[[1]])
    crick.bases <- covered[[2]]
    watson.bases <- covered[[3]]
    crick.quals <- covered[[4]]
    watson.quals <- covered[[5]]

    time <- proc.time() - ptm
    message(" ", round(time[3], 2), "s")

    ## initialize matrix with all regions as rows and covered SNVs as columns
    if (length(crick.bases) > 0 & length(watson.bases) > 0) {
        filename.IDs <- names(crick.bases)
        crickBases.m <- matrix(NA, nrow = length(filename.IDs), ncol = length(genomic.pos))
        watsonBases.m <- matrix(NA, nrow = length(filename.IDs), ncol = length(genomic.pos))
        crickQuals.m <- matrix(NA, nrow = length(filename.IDs), ncol = length(genomic.pos))
        watsonQuals.m <- matrix(NA, nrow = length(filename.IDs), ncol = length(genomic.pos))

        # message("Loading matrices")

        ## loop over all file/region IDs
        for (i in 1:length(filename.IDs)) {
            # message("Working on ", filename.IDs[i])
            filename.ID <- filename.IDs[i]
            crickBases.v <- crick.bases[[filename.ID]]
            watsonBases.v <- watson.bases[[filename.ID]]
            crickQuals.v <- crick.quals[[filename.ID]]
            watsonQuals.v <- watson.quals[[filename.ID]]

            crick.uncov.pos <- setdiff(genomic.pos, crickBases.v)
            names(crick.uncov.pos) <- rep(0, length(crick.uncov.pos))
            crickBases <- names(sort(c(crickBases.v, crick.uncov.pos)))
            crickQuals <- names(sort(c(crickQuals.v, crick.uncov.pos)))

            watson.uncov.pos <- setdiff(genomic.pos, watsonBases.v)
            names(watson.uncov.pos) <- rep(0, length(watson.uncov.pos))
            watsonBases <- names(sort(c(watsonBases.v, watson.uncov.pos)))
            watsonQuals <- names(sort(c(watsonQuals.v, watson.uncov.pos)))

            crickBases.m[i, ] <- as.numeric(crickBases)
            watsonBases.m[i, ] <- as.numeric(watsonBases)
            crickQuals.m[i, ] <- as.numeric(crickQuals)
            watsonQuals.m[i, ] <- as.numeric(watsonQuals)
        }

        ## sort matrices by the number of covered SNVs (from highest to lowest)
        cov.pos.crick <- apply(crickBases.m, 1, function(x) length(which(x > 0)))
        cov.pos.watson <- apply(watsonBases.m, 1, function(x) length(which(x > 0)))
        cov <- cov.pos.crick + cov.pos.watson

        ordered.idx <- order(cov, decreasing = T)

        if (length(ordered.idx) > 1) {
            crickBases.m <- crickBases.m[ordered.idx, ]
            watsonBases.m <- watsonBases.m[ordered.idx, ]
            crickQuals.m <- crickQuals.m[ordered.idx, ]
            watsonQuals.m <- watsonQuals.m[ordered.idx, ]
            filename.IDs <- filename.IDs[ordered.idx]
        }

        ## export matrices and file IDs as single list
        matrices <- list()
        matrices[["crick.bases"]] <- crickBases.m
        matrices[["watson.bases"]] <- watsonBases.m
        matrices[["crick.quals"]] <- crickQuals.m
        matrices[["watson.quals"]] <- watsonQuals.m
        matrices[["genomic.pos"]] <- sort(genomic.pos)
        matrices[["row.IDs"]] <- filename.IDs
        # matrices[['nonWC.reads']] <- nonWC.grl

        return(matrices)
    } else {
        matrices <- list()
        return(matrices)
    }
}
