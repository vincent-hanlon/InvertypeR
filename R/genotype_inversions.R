#' Bayesian genotyper for inversions
#'
#' Given two Strand-seq composite files (WC and WW) and a list of intervals, this computes the highest posterior probability of the possible (phased) genotypes.
#'
#' Step-by-step: It standardizes the priors so that they sum to 1, and sets prior probabilities for errors that may be present in the data. Then it chooses the binomial probabilities for
#' the Bayesian model. The function counts the forward and reverse reads in each inversion, and then genotypes them in the haploid case (details omitted: it's the same as the diploid case really).
#' The read counts and the probabilites can be combined to calculate binomial (log) likelihoods of the possible strand states of each inversion in the two composite files (e.g. WW or WC).
#' Given the strand-states we expect (accounting for the errors that may be present in the data), these likelihoods can be used to compute more (log) likelihoods: this time, for the
#' genotypes REF, HET(0|1), HET(1|0), and HOM. We convert these into regular posterior probabilties and choose the highest one, with the associated genotype.
#'
#' Note that we can phase inversions only because the WC composite file already has phased reads. This means that we know a 0|1 inversion on chr1 is on the same homolog as all 0|1 inversions
#' on chr1 in the sample, and that all chr1 1|0 inversions are on the other homolog. However, we don't know whether a 0|1 inversion on chr1 and a 0|1 inversion chr2 came from the same parent.
#' 0|1 inversions are distinguished from 1|0 inversions based on the strand switch in the WC composite file ( WC -> WW or WC -> CC).
#'
#' @param WW_reads A GRanges object (or GAlignmentPairs in the PE case) containing reads for a WW composite file. See create_composite_files() or import_bam().
#' @param WC_reads A GRanges object (or GAlignmentPairs in the PE case) containing reads for a WC composite file. See create_composite_files() or import_bam().
#' @param regions A Granges object containing genomic intervals that are thought to be inversions.
#' @param background The fraction of background reads for the WW composite file. See WWCC_background().
#' @param base_state The strand state of the WW composite file: either "WW" (mostly + reads) or "CC" (mostly - reads).
#' @param haploid_chromosomes A vector of the names of chromosomes expected to be haploid (e.g., chrX and chrY in human males). Default NULL.
#' @param prior Vector of three prior weights for inversion genotypes. For example, c("ref","het","hom") = c(0.9,0.05,0.05). Default c(0.33,0.33,0.33).
#' @param haploid_prior Vector of two prior weights for inversions on haploid chromosomes (e.g., chrX and chrY in human males). For example, c("ref", "inv") = c(0.9,0.1). Default c(0.5,0.5).
#' @return A dataframe of the regions, where each region is matched with the most probable genotype and the corresponding posterior probability, as well as some read counts.
#'
#' @export
genotype_inversions <- function(WW_reads, WC_reads, regions, background, base_state, haploid_chromosomes = NULL, prior = c(0.33, 0.33, 0.33), haploid_prior = c(0.5, 0.5)) {

    # Returning an empty dataframe if there are no regions
    if (length(regions) == 0) {
        inversions <- data.frame(character(), integer(), integer(), integer(), integer(), integer(), integer(), character(), double())
        colnames(inversions) <- c("chr", "start", "end", paste0(base_state, "_counts_W"), paste0(base_state, "_counts_C"), "WC_counts_W", "WC_counts_C", "genotype", "probability")

        return(inversions)
    }

    # Standardizing priors and splitting the heterozygous value for the two haplotypes
    prior <- c(prior[1], prior[2] / 2, prior[2] / 2, prior[3]) / sum(prior)
    # for haploid chromosomes
    haploid_prior <- haploid_prior / sum(haploid_prior)

    # Error states are (i) no error, (ii) alignment error (composite files are both WC), and (iii) deletion (nonsense strand state paires for the composite files).
    # These are effectively priors as well, and are chosen somewhat arbitrarly here.

    prior_error <- c(0.5, 0.25, 0.25)
    haploid_prior_error <- c(0.75, 0.25)

    # Have to distinguish WW and CC base state here to distinguish hom_00 and hom_01
    if (base_state == "WW") {
        theta <- c(1 - background, 0.5, 0.5, background)
    } else if (base_state == "CC") {
        theta <- c(background, 0.5, 0.5, 1 - background)
    } else {
        stop("Invalid base strand state for the WW/CC composite file. Enter WW or CC")
    }

    # read counts per region: by far the slowest step
    counts <- count_regions(WW_reads, WC_reads, regions)

    # Samples from males/haploids need to have the haploid chromosomes treated differently

    haploid_chromosomes_present <- haploid_chromosomes[haploid_chromosomes %in% as.vector(GenomicRanges::seqnames(regions))]

    if (length(haploid_chromosomes_present) > 0) {
        haploid_counts <- counts[counts[, 1] %in% haploid_chromosomes_present, ]
        counts <- counts[!counts[, 1] %in% haploid_chromosomes_present, ]
        n <- nrow(haploid_counts)

        # The binomial distribution gives us likelihoods for the read counts given the possible strand states W and C for the WWCC file
        likelihoods <- array(mapply(function(z) mapply(function(x, y) dbinom(x, size = x + y, prob = z, log = TRUE), haploid_counts[, 4], haploid_counts[, 5]), theta), dim = c(n, 4))

        # I want to give it the possibility of alignment errors (I won't allow deletions) so that it doesn't naively call 4 W reads out of 10 a 99% homozygote
        # Also dealing entirely in log-likelihoods here
        log_haploid_prior <- log(haploid_prior)
        haploid_prior_error <- log(haploid_prior_error)
        hom <- sapply(c(1:n), function(x) matrixStats::logSumExp(c(likelihoods[x, 1] + haploid_prior_error[1], likelihoods[x, 2] + haploid_prior_error[2])) + log_haploid_prior[1])
        hom_alt <- sapply(c(1:n), function(x) matrixStats::logSumExp(c(likelihoods[x, 4] + haploid_prior_error[1], likelihoods[x, 2] + haploid_prior_error[2])) + log_haploid_prior[2])
        haploid_posterior <- cbind(hom, hom_alt)

        # Bayes' Theorem
        rowsums <- sapply(c(1:n), function(x) matrixStats::logSumExp(haploid_posterior[x, ]))
        haploid_posterior <- haploid_posterior - matrix(rep(rowsums, 2), ncol = 2)

        # Converting back to regular probabilities
        haploid_posterior <- exp(haploid_posterior)
        colnames(haploid_posterior) <- c("0", "1")

        haploid_chosen <- colnames(haploid_posterior)[apply(haploid_posterior, 1, which.max)]
        haploid_chosen <- data.frame(haploid_chosen, apply(haploid_posterior, 1, max))

        # Assembling a useful dataframe to return
        haploid_inversions <- cbind(haploid_counts, haploid_chosen)
        colnames(haploid_inversions) <- c(
            "chr", "start", "end", paste0(base_state, "_counts_W"), paste0(base_state, "_counts_C"), "WC_counts_W", "WC_counts_C",
            "genotype", "probability"
        )
    } else {
        counts <- counts[!counts[, 1] %in% c("chrY"), ]
    }


    if (nrow(counts) > 0) {
        # The binomial distribution gives us likelihoods for the read counts given the possible strand states WW, WC, CC
        WW_likelihoods <- array(mapply(function(z) mapply(function(x, y) dbinom(x, size = x + y, prob = z, log = TRUE), counts[, 4], counts[, 5]), theta), dim = c(nrow(counts), 4))
        WC_likelihoods <- array(mapply(function(z) mapply(function(x, y) dbinom(x, size = x + y, prob = z, log = TRUE), counts[, 6], counts[, 7]), theta), dim = c(nrow(counts), 4))
        n <- nrow(WW_likelihoods)

        # computing the likelihoods for each inversion genotype, under tha assumption that it can either appear plainly, with alignment errors, or with an additional heterozgous deletion
        # The prior term for the deletion needs to be divided by two, since it can be calculated two ways
        # Everything written in log-arithmetic to avoid floating point issues.
        prior_error <- log(prior_error)
        log_prior <- log(prior)
        hom <- sapply(c(1:n), function(x) {
            matrixStats::logSumExp(c(WW_likelihoods[x, 1] + WC_likelihoods[x, 1] - log(2) + prior_error[3], WW_likelihoods[x, 1] + WC_likelihoods[x, 4] - log(2) +
                prior_error[3], WC_likelihoods[x, 2] + prior_error[1] + WW_likelihoods[x, 1], WW_likelihoods[x, 2] + WC_likelihoods[x, 2] + prior_error[2])) + log_prior[1]
        })

        hom_alt <- sapply(c(1:n), function(x) {
            matrixStats::logSumExp(c(WW_likelihoods[x, 4] + WC_likelihoods[x, 1] - log(2) + prior_error[3], WW_likelihoods[x, 4] + WC_likelihoods[x, 4]
                - log(2) + prior_error[3], WC_likelihoods[x, 2] + prior_error[1] + WW_likelihoods[x, 4], WW_likelihoods[x, 2] + WC_likelihoods[x, 2] + prior_error[2])) + log_prior[4]
        })

        het_01 <- sapply(c(1:n), function(x) {
            matrixStats::logSumExp(c(WW_likelihoods[x, 2] + WC_likelihoods[x, 1] + prior_error[1], WW_likelihoods[x, 2] + WC_likelihoods[x, 2] +
                prior_error[2], WC_likelihoods[x, 1] + WW_likelihoods[x, 1] - log(2) + prior_error[3], WC_likelihoods[x, 1] + WW_likelihoods[x, 4] - log(2) + prior_error[3])) +
                log_prior[2]
        })

        het_10 <- sapply(c(1:n), function(x) {
            matrixStats::logSumExp(c(WW_likelihoods[x, 2] + WC_likelihoods[x, 4] + prior_error[1], WW_likelihoods[x, 2] + WC_likelihoods[x, 2] +
                prior_error[2], WC_likelihoods[x, 4] + WW_likelihoods[x, 1] - log(2) + prior_error[3], WC_likelihoods[x, 4] + WW_likelihoods[x, 4] - log(2) + prior_error[3])) +
                log_prior[3]
        })

        # Computing the posterior probabilities for the four phased genotypes based on Bayes Theorem
        posterior <- cbind(hom, het_01, het_10, hom_alt)
        rowsums <- sapply(c(1:n), function(x) matrixStats::logSumExp(posterior[x, ]))
        posterior <- posterior - matrix(rep(rowsums, 4), ncol = 4)

        # Converting back to regular probabilities
        posterior <- exp(posterior)
        colnames(posterior) <- c("0|0", "0|1", "1|0", "1|1")

        # Then we select the genotype (and posterior probability) that is most likely
        chosen <- colnames(posterior)[apply(posterior, 1, which.max)]
        chosen <- data.frame(chosen, apply(posterior, 1, max))

        # Assembling a useful dataframe to return
        inversions <- cbind(counts, chosen)
        colnames(inversions) <- c("chr", "start", "end", paste0(base_state, "_counts_W"), paste0(base_state, "_counts_C"), "WC_counts_W", "WC_counts_C", "genotype", "probability")
    }

    # Combining the results from the sex/haploid chromosomes and the autosomes/diploid chromosomes appropriately (in the haploid case)
    if (length(haploid_chromosomes_present) > 0 & nrow(counts) > 0) {
        inversions <- rbind(inversions, haploid_inversions)
    } else if (nrow(counts) <= 0 & length(haploid_chromosomes_present) > 0) {
        inversions <- haploid_inversions
    }

    # sorting the inversions by posterior probability
    inversions <- inversions[order(-inversions$probability), ]

    return(inversions)
}
