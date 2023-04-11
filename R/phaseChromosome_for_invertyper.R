#' A stripped-down version of the StrandPhaseR function phaseChromosome by David Porubsky
#'
#' This is a simplified function that removes some of the functionality (and runtime) of David's original StrandPhaseR function
#' The arguments are the same (or effectively the same) as David's corresponding function.
#'
#' The descriptions of arguments here are just copied from StrandPhaseR
#'
#' @param input_folder Path to the bam files to process
#' @param output_folder Output directory.
#' @param positions Filename with listed position of SNVs for given chromosome (format: chrName SNVpos).
#' @param WCregions Filename of all WC region for a given chromosome (format: chrName:Start:End:FileName).
#' @param chromosome If only a subset of the chromosomes should be processed, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param min.baseq Minimum base quality to consider a base for phasing.
#' @param num.iterations Number of iteration to sort watson and crick matrices.
#' @param translateBases ...
#' @param fillGaps ...
#' @param splitPhasedReads Set to \code{TRUE} if you want to split reads per haplotype.
#' @param compareSingleCells Set to \code{TRUE} if you want to compare haplotypes at the single-cell level.
#' @param exportVCF ...
#'
#' @return A subset of the StrandPhaseR output that phases the strand state of wc regions
#' 
#' @export
phaseChromosome_for_invertyper <- function(input_folder, output_folder='./StrandPhaseR_analysis', positions=NULL, WCregions=NULL, chromosome=NULL, pairedEndReads=TRUE, min.mapq=10, min.baseq=30, num.iterations=2, translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, compareSingleCells=FALSE, exportVCF=NULL) {

        GenomeInfoDb::seqlevels(positions, pruning.mode='coarse') <- chromosome
        GenomeInfoDb::seqlevels(WCregions, pruning.mode='coarse') <- chromosome
  #load data into matrix
  matrices <- loadMatrices_for_invertyper(input_folder=input_folder, positions=positions, WCregions=WCregions, pairedEndReads=pairedEndReads, min.mapq=min.mapq, min.baseq=min.baseq)
  #Check if sufficient data were loaded
  if (length(matrices) > 0) {
    #phase data
    srt.matrices <-  StrandPhaseR::sortMatrices(data.object=matrices, num.iterations=num.iterations)
    assem.haps <-  StrandPhaseR::assembleHaps(data.object=srt.matrices, translateBases=translateBases)
    hap1 <- data.frame(names(assem.haps$hap1.files), do.call(rbind, lapply(assem.haps$hap1.files, rbind)))
    names(hap1) <- c("Filenames", "Simil", "Disimil")

    return(hap1$Filenames)
 }
} 



