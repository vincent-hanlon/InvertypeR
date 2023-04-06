#' A modified implementation of the core StrandPhaseR function, "strandPhaseR"
#' 
#' This function is effectively a stripped-down version of David Porubsky's strandPhaseR, but with some additional post-processing of my own to return a nice GRanges of phased WC regions. 
#' It sneakily uses the galignmentlist (BAM files loaded into R) rather then reading them from the files all over again. This is accomplished by overwriting some StrandPhaseR functions.
#'
#' @param numCPU An integer, the number of threads to use
#' @param positions The path to a VCF file containing heterozygous SNPs 
#' @param WCregions A pre-processed GRanges object that contains intervals in which specified libraries have strand-state WC (Watson-Crick)
#' @param chromosomes A character vector of chromosome names to be used
#' @param paired_reads Boolean. Are the reads paired-end?
#' @param num.iterations An integer. See StrandPhaseR for details. The default here is usually fine.
#' @param galignmentslist A list of GenomicAlignments objects containing the input Strand-seq libraries
#'
#' @return A GRanges object such that the input WCregions are annotated with 'wc' or 'cw' according to their phase
#'
#' @export
strandPhaseR_for_invertyper <- function(numCPU=4, positions=NULL, WCregions=NULL, chromosomes=NULL, paired_reads=TRUE, num.iterations=3, galignmentslist=galignmentslist) {

`%dopar%` <- foreach::`%dopar%`
	
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

	output_folder<-'./SPR_output'
	input_folder<-'./'

	if(is.null(conf[['chromosomes']]) && !is.null(conf[['WCregions']])){
		conf[['chromosomes']] <- sort(as.character.factor(unique(GenomicRanges::seqnames(conf[['WCregions']]))))
	}

	## Loading in list of SNV positions and locations of WC regions
	snvs <- suppressWarnings(suppressMessages(StrandPhaseR::vcf2ranges(vcfFile=conf[['positions']], genotypeField=1, chromosome=conf[['chromosomes']])))

	conf[['chromosomes']] <- as.character(conf[['chromosomes']])


	all_phased_WCregions <- list()
	
	# need to awkwardly make this accessible to all sub-functions...
	galignmentslist_global_for_invertyper <<- galignmentslist


message("bamregions should write a whole bunch of empty files if the correct one is being used")
	if(conf[['numCPU']]>1) {


    cl <- parallel::makeCluster(conf[['numCPU']])
    doParallel::registerDoParallel(cl)
parallel::clusterExport(cl=cl, 'galignmentslist_global_for_invertyper')  
    all_phased_WCregions <- foreach::foreach (i = conf[['chromosomes']], .packages=c('StrandPhaseR', 'invertyper')) %dopar% {
    
        tC <- tryCatch({

message("remember to add suppressMessages( back to strandphaser here!") 
message( "using the right foreach loop")
phaseChromosome_for_invertyper(input_folder=input_folder, 
output_folder=output_folder, positions=snvs[GenomicRanges::seqnames(snvs) == i],
                                        chromosome=i,
                                       WCregions=WCregions[GenomicRanges::seqnames(WCregions)==i], pairedEndReads=conf[['pairedEndReads']], min.mapq=conf[['min.mapq']],
                                min.baseq=20, num.iterations=conf[['num.iterations']],
                                       translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, compareSingleCells=FALSE, exportVCF=NULL)


	 }, error = function(err) {
          stop(i,'\n',err)
        })
	}		
    parallel::stopCluster(cl)


	} else {


	for(i in conf[['chromosomes']]){

	message("remember to add suppressMessages( back to strandphaser!")
	temp <- phaseChromosome_for_invertyper(input_folder=input_folder, output_folder=output_folder, positions=snvs[GenomicRanges::seqnames(snvs) == i], 
					chromosome=i,
                                       WCregions=WCregions[GenomicRanges::seqnames(WCregions)==i], pairedEndReads=conf[['pairedEndReads']], min.mapq=conf[['min.mapq']], 
				min.baseq=20, num.iterations=conf[['num.iterations']],
                                       translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, compareSingleCells=FALSE, exportVCF=NULL)
	all_phased_WCregions <- append(all_phased_WCregions, temp)
	}

	}



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




