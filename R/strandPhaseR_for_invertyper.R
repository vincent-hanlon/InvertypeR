
#' A modified implementation of the core StrandPhaseR function, "strandPhaseR"
#' 
#' This function is effectively a stripped-down version of David Porubsky's strandPhaseR, but with some additional post-processing of my own to return a nice GRanges of phased WC regions. 
#' It sneakily uses the galignmentlist (BAM files loaded into R) rather then reading them from the files all over again. This is accomplished by overwriting some StrandPhaseR functions.
#'
#' @export
strandPhaseR_for_invertyper <- function(numCPU=4, positions=NULL, WCregions=NULL, chromosomes=NULL, paired_reads=TRUE, num.iterations=3, galignmentslist=galignmentslist) {

R.utils::reassignInPackage("bamregion2GRanges", "StrandPhaseR", bamregion2GRanges_for_invertyper)
	
	### Helper functions ###
	#=======================
	as.object <- function(x) {
	  return(eval(parse(text=x)))
	}

        extract_wc <- function(x){

                x <- unlist(strsplit(x,"__"))
                x[5] <- x[3]
                y <- unlist(strsplit(x[2],":|-"))
                x[2] <- y[1]
                x[3] <- y[2]
                x[4] <- y[3]

                return(x)
        }


	## Put options into list and merge with conf
	params <- list(numCPU=numCPU, positions=positions, WCregions=WCregions, chromosomes=chromosomes, pairedEndReads=paired_reads, num.iterations=num.iterations)
	conf <- params[setdiff(names(params),names(NULL))]

	#===================
	### Input checks ###
	#===================
	## Check user input

	outputfolder<-'./SPR_output'
	inputfolder<-'./'

	if(is.null(conf[['chromosomes']]) && !is.null(conf[['WCregions']])){
		conf[['chromosomes']] <- sort(as.character.factor(unique(GenomicRanges::seqnames(conf[['WCregions']]))))
	}

	## Loading in list of SNV positions and locations of WC regions
	snvs <- suppressWarnings(suppressMessages(StrandPhaseR::vcf2ranges(vcfFile=conf[['positions']], genotypeField=1, chromosome=conf[['chromosomes']])))

	conf[['chromosomes']] <- as.character(conf[['chromosomes']])


	all_phased_WCregions <- list()
	
	# need to awkwardly make this accessible to all sub-functions...
	galignmentslist_global_for_invertyper <<- galignmentslist

	for(i in conf[['chromosomes']]){

message("remember to add suppressMessages( back to strandphaser!")
	temp <- phaseChromosome_for_invertyper(inputfolder=inputfolder, outputfolder=outputfolder, positions=snvs[GenomicRanges::seqnames(snvs) == i], 
					chromosome=i,
                                       WCregions=WCregions[GenomicRanges::seqnames(WCregions)==i], pairedEndReads=conf[['pairedEndReads']], min.mapq=conf[['min.mapq']], min.baseq=20, num.iterations=conf[['num.iterations']],
                                       translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, compareSingleCells=FALSE, exportVCF=NULL)
	all_phased_WCregions <- append(all_phased_WCregions, temp)
	}

#	all_phased_WCregions <- suppressMessages(lapply(conf[['chromosomes']], StrandPhaseR::phaseChromosome, inputfolder=inputfolder, outputfolder=outputfolder, positions=snvs[GenomicRanges::seqnames(snvs) == conf[['chromosomes']]], 
#					WCregions=WCregions, pairedEndReads=conf[['pairedEndReads']], min.mapq=conf[['min.mapq']], min.baseq=20, num.iterations=conf[['num.iterations']], 
#					translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, compareSingleCells=FALSE, exportVCF=NULL))#, mc.cores=numCPU))

	rm('galignmentslist_global_for_invertyper', envir=globalenv())
	invisible(gc())

	all_phased_WCregions <- lapply(unlist(all_phased_WCregions),extract_wc)

        all_phased_WCregions <- as.data.frame(do.call(rbind,all_phased_WCregions))
        all_phased_WCregions[all_phased_WCregions[5]=="W",5] <- "wc"
        all_phased_WCregions[all_phased_WCregions[5]=="C",5] <- "cw"

        names(all_phased_WCregions) <- c("filename","chr","start","end","states")

        all_phased_WCregions <- GenomicRanges::makeGRangesFromDataFrame(all_phased_WCregions[,c(2,3,4,1,5)], keep.extra.columns=TRUE)
        GenomeInfoDb::seqlevels(all_phased_WCregions) <- GenomeInfoDb::seqlevels(galignmentslist[[1]])

        return(all_phased_WCregions)

}




