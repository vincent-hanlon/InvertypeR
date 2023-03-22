#' A stripped-down version of the StrandPhaseR function phaseChromosome
#'
#' This is a simplified function that removes some of the functionality (and runtime) of David Porubsky's original StrandPhaseR function
#'
#' @export
phaseChromosome_for_invertyper <- function(inputfolder, outputfolder='./StrandPhaseR_analysis', positions=NULL, WCregions=NULL, chromosome=NULL, pairedEndReads=TRUE, min.mapq=10, min.baseq=30, num.iterations=2, translateBases=TRUE, fillMissAllele=NULL, 
splitPhasedReads=FALSE, compareSingleCells=FALSE, exportVCF=NULL) {

	message("Using Vincents altered PhaseChromosome...")
        GenomeInfoDb::seqlevels(positions, pruning.mode='coarse') <- chromosome
        GenomeInfoDb::seqlevels(WCregions, pruning.mode='coarse') <- chromosome
  #load data into matrix
  matrices <- StrandPhaseR::loadMatrices(inputfolder=inputfolder, positions=positions, WCregions=WCregions, pairedEndReads=pairedEndReads, min.mapq=min.mapq, min.baseq=min.baseq)

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



