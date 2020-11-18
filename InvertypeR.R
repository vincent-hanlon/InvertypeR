#Next up: InvertypeR as an R package
#This set of functions implements a binomial Bayesian model to genotype putative inversions
#It can also adjust the endpoints of called inversions based on deltaW values in the composite files





library(GenomicAlignments)

###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Wrapper functon for the other functions below.
#Note that I am assuming 1-based coordinates

#WW_bam, WC_bam are the paths to the two composite BAM files
#bed is a bed file containing the coordinates of the putative inversions
#blacklist is a bed file containing the coordinates of regions with unreliable strand-seq data
#prior and prior_male are vectors of prior weights defined by the user. prior=c(P(ref_hom), P(het), P(alt_hom)) and prior_male=c(P(ref), P(alt))
#adjust_method should be one of "raw" (no adjustment of endpoints for inversions), "merge" (merge inversions of the same genotype), "deltas" (adjust based on deltaW values if 
	#possible, or if not merge, or if not take the largest of a set of overlapping intervals), "minimal" (for each set of overlapping inversions of the same genotype, take the interval common to 
	#all of them if one exists), or "all" (does "deltas" and then also attempts to find strand switches near to genotypes with posterior < confidence). Note that the adjust method requires 
	#that the new inteval overlaps at least one of the overlapping inversions it came from, so inversions can't "walk" away from their original coordinates. If after adjustment 
	#new overlaps are formed (i.e.adjacent intervals expands so they overlap), then they are merged or the largest is taken, as above.

invertyper <- function(WW_bam, WC_bam, bed, blacklist="", paired_reads=TRUE, sex="female", confidence=0.95,prior=c(0.333,0.333,0.333), prior_male=c(0.5,0.5), 
					output_file="inversions.txt", adjust_method=c("raw","merge","deltas","minimal", "all", "low")){

	#Loading the BED files as GRanges objects
	regions <- import(bed)

	#Takes reads from the BAM files that overlap the regions and stores them as GRanges
	reads <- read_regions(WW_bam, WC_bam, regions, paired_reads=paired_reads, blacklist=blacklist)

        #Accurate background estimate, plus base strand state for the WW/CC file
        base <- WWCC_background(WW_bam, binsize=1000000, paired_reads=paired_reads)

	#Genotyping the inversions
	inversions <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=regions, background=base[[1]], base_state=base[[2]],  sex=sex, 
		prior=prior, prior_male=prior_male)


	#two methods under development for dealing with overlapping intervals
	if( adjust_method != "raw" ){
		
		inversions <- adjust_intervals(inversions=inversions, reads=reads, confidence=confidence, base=base,  sex=sex, prior=prior, prior_male=prior_male, 
			adjust_method=adjust_method, paired_reads=paired_reads,blacklist=blacklist)
	}

	#Marking inversions with <0.001 reads per bp as being low-density, because they may span centromeres or chromosome ends/
	#The 0.001 figure is somewhat arbitrary, but it's adjusted based on the read count in the WW file of any submitted inversions. 
	inversions[inversions[,9]>=confidence & inversions[,8]!=0 & inversions[,8]!="0|0" & rowSums(inversions[,c(4:7)])/(inversions[,3]-inversions[,2] + 1) < (0.001*base[[3]]/27734969), 10] <- 
			"low read density: check for inaccurate inversion coordinates"

	colnames(inversions)[10] <- "low_read_density"

	write.table(inversions, output_file, quote=FALSE, sep="\t", row.names=FALSE)

	return(inversions)
}

##########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#A function to get a more precise estimate of the background associated with a WW or CC composite BAM file
#Also checks whether the ground state is WW or CC (relative to the definitions in the strand_states fn) and returns an error if it looks like WC/CW

WWCC_background <- function(WW_bam, binsize=1000000, paired_reads=TRUE){

	file <- BamFile(WW_bam)
	#chromosome lengths for the first 22 chromosomes 
	chr_lengths <- scanBamHeader(file)$targets[1:22]
	#Creating genomic bins, default size 1 Mb
	bins <-	tileGenome(chr_lengths, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)

	#Generating parameters for countBam that are appropriate for pe or se reads
	if(paired_reads) {

		p_minus <- ScanBamParam(flag = scanBamFlag(isProperPair=TRUE,isUnmappedQuery=FALSE,isDuplicate=FALSE,isFirstMateRead=TRUE, isMinusStrand=TRUE),mapqFilter=10, which=bins)
		p_plus <- ScanBamParam(flag = scanBamFlag(isProperPair=TRUE,isUnmappedQuery=FALSE,isDuplicate=FALSE,isFirstMateRead=TRUE, isMinusStrand=FALSE),mapqFilter=10, which=bins)	
		
        } else {

		p_minus <- ScanBamParam(flag = scanBamFlag(isPaired=FALSE,isUnmappedQuery=FALSE,isDuplicate=FALSE, isMinusStrand=TRUE),mapqFilter=10, which=bins)
		p_plus <- ScanBamParam(flag = scanBamFlag(isPaired=FALSE,isUnmappedQuery=FALSE,isDuplicate=FALSE, isMinusStrand=FALSE),mapqFilter=10, which=bins)

        }

	#Counting reads by strand per bin, and using this to estimate background per bin
	num_c_bins <- countBam(file, param=p_minus)[,"records"]
	num_w_bins <- countBam(file, param=p_plus)[,"records"]

	num_reads <- sum(num_c_bins+num_w_bins)

	background_bins <- num_c_bins/(num_c_bins + num_w_bins)
	background_bins <- background_bins[!is.na(background_bins)]
	
	#Kernel density of the background_bins distribution 
	background_density <- density(background_bins)
	#Finding the peak of the distribution and reporting that as background, to avoid always-WC regions, inversions, etc.
	i <- which.max(background_density$y)	
	background <- background_density$x[i]

	
	#Figuring out whether we have WW or CC so the base state can be matched with the strand_state() output
	if( background < 0.1 ) {

               base_state <- "WW"

        } else if( background > 0.9 ) {

               base_state <- "CC"

        } else if( ( background > 0.3) & (background < 0.7) ) {

                stop("the input file seems to be WC/CW, not WW/CC!")

        } else {

                stop("the input file has >10% background!")

        }


	return(list(min(background, 1-background),base_state,num_reads))

}

###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Jointly estimates the genotype of a region based on the WC and WW composite files.
#Note that the phasing is only relative to other inversions on the same chromosome, and 1st and 2nd homolog are arbitrarily identified according to the base state (WW or CC)


genotype_inversions <- function(WW_reads, WC_reads, regions, background, base_state, sex="female", prior=c(0.33,0.33,0.33), prior_male=c(0.5,0.5)){

	#Returning an empty dataframe if there are no regions
	if(length(regions)==0){ 

		inversions <- data.frame(character(),integer(),integer(),integer(),integer(),integer(),integer(),character(),double())
                colnames(inversions) <- c("chr","start","end",paste0(base_state,"_counts_W"), paste0(base_state,"_counts_C"), "WC_counts_W", "WC_counts_C","genotype","probability")

		return(inversions)
	 }

	#Standardizing priors and splitting the heterozygous value for the two haplotypes
	prior <- c(prior[1],prior[2]/2,prior[2]/2, prior[3])/sum(prior)
	#for the sex chromosomes in males
	prior_male <- prior_male/sum(prior_male)

	#Error states are (i) no error, (ii) alignment error (composite files are both WC), and (iii) deletion (nonsense strand state paires for the composite files).
	#These are effectively priors as well, and are chosen somewhat arbitrarly here.  

        prior_error <- c(0.5, 0.25, 0.25)
 	prior_error_male <- c(0.75, 0.25)


	#Have to distinguish WW and CC base state here to distinguish hom_00 and hom_01
        if(base_state=="WW"){

                theta <- c(1-background,0.5,0.5,background)

        } else if(base_state=="CC") {

                theta <- c(background,0.5,0.5,1-background)

        } else {

                stop("Invalid base strand state for the WW/CC composite file. Enter WW or CC")
        }

	
	#read counts per region: by far the slowest step
	counts <- count_regions(WW_reads, WC_reads, regions)	

	sex_chromosomes <- NULL

	#Samples from males need to have the haploid sex chromosomes treated differently
	if(sex == "male"){

		#Figuring out which sex chromosomes are present
		sex_chromosomes <- intersect(unique(as.character(seqnames(regions))), c("chrX", "chrY"))

		if(length(sex_chromosomes) > 0) {

			male_counts <- counts[counts[,1] %in% sex_chromosomes,]
			counts <- counts[!counts[,1] %in% sex_chromosomes,]
			n <- nrow(male_counts)			

		        #The binomial distribution gives us likelihoods for the read counts given the possible strand states W and C for the WWCC file
        		likelihoods <- array(mapply(function(z) mapply(function(x,y) dbinom(x,size=x+y,prob=z,log=TRUE),male_counts[,4],male_counts[,5]), theta),
				dim=c(n,4))

			#I want to give it the possibility of alignment errors (I won't allow deletions) so that it doesn't naively call 4 W reads out of 10 a 99% homozygote		
			#Also dealing entirely in log-likelihoods here 	
			log_prior_male <- log(prior_male)
			prior_error_male <- log(prior_error_male)
			hom <- sapply(c(1:n), function(x) logSumExp(c(likelihoods[x,1] + prior_error_male[1], likelihoods[x,2] + prior_error_male[2])) + log_prior_male[1])
        	        hom_alt <- sapply(c(1:n), function(x) logSumExp(c(likelihoods[x,4] + prior_error_male[1], likelihoods[x,2] + prior_error_male[2])) + log_prior_male[2])
			male_posterior <- cbind(hom,hom_alt)

			#Bayes' Theorem
			rowsums <- sapply(c(1:n), function(x) logSumExp(male_posterior[x,]))
                	male_posterior <- male_posterior - matrix(rep(rowsums, 2), ncol=2)

                	#Converting back to regular probabilities
               		male_posterior <- exp(male_posterior)
		        colnames(male_posterior) <- c("0","1")
			
			male_chosen <- colnames(male_posterior)[apply(male_posterior,1,which.max)]
	        	male_chosen <- data.frame(male_chosen, apply(male_posterior,1,max))

	        	#Assembling a useful dataframe to return
        		male_inversions <- cbind(male_counts, male_chosen)
			colnames(male_inversions) <- c("chr","start","end",paste0(base_state,"_counts_W"), paste0(base_state,"_counts_C"), 
				"WC_counts_W", "WC_counts_C","genotype","probability")
		}	

	} else {
	
		sex_chromosomes <- "chrY"
		counts <- counts[!counts[,1] %in% sex_chromosomes,]

	}


	if(nrow(counts)>0) {


		#The binomial distribution gives us likelihoods for the read counts given the possible strand states WW, WC, WW	        
		WW_likelihoods <- array(mapply(function(z) mapply(function(x,y) dbinom(x,size=x+y,prob=z,log=TRUE),counts[,4],counts[,5]), theta),dim=c(nrow(counts),4))
		WC_likelihoods <- array(mapply(function(z) mapply(function(x,y) dbinom(x,size=x+y,prob=z,log=TRUE),counts[,6],counts[,7]), theta),dim=c(nrow(counts),4))
		n <- nrow(WW_likelihoods)

		#computing the likelihoods for each inversion genotype, under tha assumption that it can either appear plainly, with alignment errors, or with an additional heterozgous deletion
		#The prior term for the deletion needs to be divided by two, since it can be calculated two ways
		#Everything written in log-arithmetic to avoid floating point issues.
		prior_error <- log(prior_error)
		log_prior <- log(prior)
		hom <- sapply(c(1:n), function(x) logSumExp(c(WW_likelihoods[x,1] + WC_likelihoods[x,1] -log(2) + prior_error[3], WW_likelihoods[x,1] + WC_likelihoods[x,4] -log(2) + 
			prior_error[3], WC_likelihoods[x,2] + prior_error[1] + WW_likelihoods[x,1], WW_likelihoods[x,2] + WC_likelihoods[x,2] + prior_error[2])) + log_prior[1])

		hom_alt <-  sapply(c(1:n), function(x) logSumExp(c(WW_likelihoods[x,4] + WC_likelihoods[x,1] -log(2) + prior_error[3], WW_likelihoods[x,4] + WC_likelihoods[x,4] 
			-log(2) + prior_error[3], WC_likelihoods[x,2] + prior_error[1] + WW_likelihoods[x,4], WW_likelihoods[x,2] + WC_likelihoods[x,2] + prior_error[2])) + log_prior[4])

		het_01 <- sapply(c(1:n), function(x) logSumExp(c(WW_likelihoods[x,2] + WC_likelihoods[x,1] + prior_error[1], WW_likelihoods[x,2] + WC_likelihoods[x,2] + 
			prior_error[2], WC_likelihoods[x,1] + WW_likelihoods[x,1] - log(2) + prior_error[3],  WC_likelihoods[x,1] + WW_likelihoods[x,4] - log(2) + prior_error[3])) + 
			log_prior[2])

		het_10 <- sapply(c(1:n), function(x) logSumExp(c(WW_likelihoods[x,2] + WC_likelihoods[x,4] + prior_error[1], WW_likelihoods[x,2] + WC_likelihoods[x,2] + 
			prior_error[2], WC_likelihoods[x,4] + WW_likelihoods[x,1] - log(2) + prior_error[3],  WC_likelihoods[x,4] + WW_likelihoods[x,4] - log(2) + prior_error[3])) + 
			log_prior[3])

		#Computing the posterior probabilities for the four phased genotypes based on Bayes Theorem
		posterior <- cbind(hom, het_01, het_10, hom_alt)
		rowsums <- sapply(c(1:n), function(x) logSumExp(posterior[x,]))
		posterior <- posterior - matrix(rep(rowsums, 4), ncol=4)
		
		#Converting back to regular probabilities
		posterior <- exp(posterior)
	        colnames(posterior) <- c("0|0","0|1","1|0","1|1")

		#Then we select the genotype (and posterior probability) that is most likely
	        chosen <- colnames(posterior)[apply(posterior,1,which.max)]
	        chosen <- data.frame(chosen, apply(posterior,1,max))

       		#Assembling a useful dataframe to return
      		inversions <- cbind(counts, chosen)
      	 	colnames(inversions) <- c("chr","start","end",paste0(base_state,"_counts_W"), paste0(base_state,"_counts_C"), "WC_counts_W", "WC_counts_C","genotype","probability")


	}




	#Combining the results from the sex chromosomes and the autosomes appropriately (in the male case)	
	if(sex == "male" & length(sex_chromosomes) > 0 & nrow(counts) > 0) {
		
		inversions <- rbind(inversions, male_inversions)
	
	} else if(nrow(counts) <= 0 & sex == "male" & length(sex_chromosomes) > 0) {
	
		inversions <- male_inversions
	}


	#sorting the inversions by posterior probability
	inversions <- inversions[order(-inversions$probability),]

        return(inversions)

}

########################################################################################################################################################################################

#By far the slowest part of the genotype_intervals function
#vectorizing using lapply didn't speed this up

count_regions <- function(WW_reads, WC_reads, regions) {

       #Extracting counts of W and C reads by region
        WW_w <- countOverlaps(subject=WW_reads[strand(WW_reads)=="+"],query=regions)
        WW_c <- countOverlaps(subject=WW_reads[strand(WW_reads)=="-"],query=regions)
        WC_w <- countOverlaps(subject=WC_reads[strand(WC_reads)=="+"],query=regions)
        WC_c <- countOverlaps(subject=WC_reads[strand(WC_reads)=="-"],query=regions)

        #Assemmbling a data frame with region info
        counts <- data.frame(seqnames(regions), start(regions), end(regions), WW_w, WW_c, WC_w, WC_c)
        colnames(counts) <- c("chr","start","end", "WW_w", "WW_c", "WC_w","WC_c")

        return(counts)
}

###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Reads in intervals within strandseq composite files
#Didn't want to have to duplicate this bit of code in the genotype_intervals() function

read_regions <- function(WW_bam, WC_bam, regions, paired_reads=TRUE, blacklist=blacklist) {

	if(blacklist!=""){
        	regions <- setdiff(reduce((regions+1000000),min.gapwidth=1000),import(blacklist))
        } else {
        	regions <- reduce((regions+1000000),min.gapwidth=1000)
        }


        if(paired_reads) {

		#Reading in files. Warnings suppressed because readGAlignment... generates them for no good reason
		#Note that the intervals are extended by 1Mb (so that we can read in flanking regions for other purposes) and reduced so they aren't overlapping
                p1 <- ScanBamParam(flag = scanBamFlag(isPaired=TRUE,isUnmappedQuery=FALSE,isDuplicate=FALSE),
			mapqFilter=10,which=regions,what=c("mapq"))
                WW_reads <- suppressWarnings(readGAlignmentPairs(WW_bam,param=p1))
                WC_reads <- suppressWarnings(readGAlignmentPairs(WC_bam,param=p1))

        } else {

                p1 <- ScanBamParam(flag = scanBamFlag(isPaired=FALSE,isUnmappedQuery=FALSE,isDuplicate=FALSE),mapqFilter=10,which=regions,what=c("mapq"))
                WW_reads <- suppressWarnings(readGAlignments(WW_bam,param=p1))
                WC_reads <- suppressWarnings(readGAlignments(WC_bam,param=p1))

         }

return(list(WW_reads,WC_reads))

}
 
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Code for pairing two fragments of a PE GRanges object taken directly from breakpointR for convenience. I did not write this code!
#Used so that deltaWCalculator can be applied for PE reads too

pair2frgm = function(data.raw){

        data.prop.pairs <- data.raw[GenomicAlignments::isProperPair(data.raw)]
        data.first <- as(GenomicAlignments::first(data.prop.pairs), 'GRanges')
        data.last <- as(GenomicAlignments::last(data.prop.pairs), 'GRanges')

        data.first.plus <- data.first[strand(data.first) == '+']
        data.first.minus <- data.first[strand(data.first) == '-']
        data.last.plus <- data.last[strand(data.last) == '+']
        data.last.minus <- data.last[strand(data.last) == '-']

        frag.plus.mapq <- data.first.plus$mapq + data.last.minus$mapq
        frag.minus.mapq <- data.first.minus$mapq + data.last.plus$mapq

        data.frag.plus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.plus), ranges=IRanges(start=start(data.first.plus), end=end(data.last.minus)), strand=strand(data.first.plus), mapq=frag.plus.mapq)
        GenomeInfoDb::seqlengths(data.frag.plus) <- GenomeInfoDb::seqlengths(data.first)
        data.frag.minus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.minus), ranges=IRanges(start=start(data.last.plus), end=end(data.first.minus)), strand=strand(data.first.minus), mapq=frag.minus.mapq)
        GenomeInfoDb::seqlengths(data.frag.minus) <- GenomeInfoDb::seqlengths(data.first)


        data <- GenomicRanges::sort(c(data.frag.plus, data.frag.minus), ignore.strand=TRUE)


        return(data)
}

############################################################################################################################################################################################

#Processes reads for deltaW adjustment

process_reads <- function(WW_reads, WC_reads, paired_reads=TRUE){

        #Processing input reads for Aaron's deltaWCalculator

        if (paired_reads) {

                WW_reads <- pair2frgm(WW_reads)
                WC_reads <- pair2frgm(WC_reads)


        } else {

                WW_reads <- GRanges(seqnames=seqnames(WW_reads),ranges=IRanges(start=start(WW_reads), end=end(WW_reads)), strand=strand(WW_reads), mcols(WW_reads), seqlengths=seqlengths(WW_reads))
                WC_reads <- GRanges(seqnames=seqnames(WC_reads),ranges=IRanges(start=start(WC_reads), end=end(WC_reads)), strand=strand(WC_reads), mcols(WC_reads), seqlengths=seqlengths(WC_reads))

        }

        #Calculating deltaW values
        WW_d <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 5))
        WC_d <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 5))

        WW_d_10 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 10))
        WC_d_10 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 10))

        WW_d_20 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 20))
        WC_d_20 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 20))

        WW_d_40 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 40))
        WC_d_40 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 40))

        WW_d_80 <- suppressWarnings(breakpointR::deltaWCalculator(WW_reads, 80))
        WC_d_80 <- suppressWarnings(breakpointR::deltaWCalculator(WC_reads, 80))

        #Taking an arbitrary sum of deltaWs to look near and far (including 2x the deltaW10)
        WW_d$deltaW <- WW_d$deltaW + 2*WW_d_10$deltaW + WW_d_20$deltaW + WW_d_40$deltaW + WW_d_80$deltaW
        WC_d$deltaW <- WC_d$deltaW + 2*WC_d_10$deltaW + WC_d_20$deltaW + WC_d_40$deltaW + WC_d_80$deltaW



        return(list(WW_reads, WC_reads,WW_d,WC_d))

}



############################################################################################################################################################################################

#Adjusts the breakpoints of inversions

adjust_deltaW <- function(WW_reads, WC_reads, WW_d, WC_d, inversions, genotype=c("het","hom","low")){


        #Next I want to merge and expand the inversions. The wider, expanded inversion intervals will be used to look for alternative breakpoints
        inversions <- sort(reduce(inversions))
        #It's simpler to get rid of strands----otherwise, precede and follow misbehave
        strand(WW_reads) <- "*"

        #the width of inversions and hom in reads
        num_inversions <- pmax(countOverlaps(inversions, WW_reads),10)

        #Follow and precede don't work very well, especially not with empty GRanges, so I need to force zeros in the empty case
        #And also add a few annoying steps to deal with NAs appropriately.
        r <- data.frame(WW_reads)

        if(length(inversions) == 0){

                f_inversions <- 0
                p_inversions <- 0

        } else {

                f_inversions <- follow(inversions,WW_reads)
                p_inversions <- precede(inversions,WW_reads)

		if(length(inversions[is.na(f_inversions)]) > 0){
			f_inversions[is.na(f_inversions)] <- precede(inversions[is.na(f_inversions)], WW_reads)
		}
		
                if(length(inversions[is.na(p_inversions)]) > 0){
                        p_inversions[is.na(p_inversions)] <- follow(inversions[is.na(p_inversions)], WW_reads)
                }
	

        }

        #Calculating the start and end coordinates for the widened interval, and making sure they don't spill over onto the next chromosome
        #The new intervals are usually 3x as wide as the originals
        start_inversions <- ifelse(as.character(seqnames(WW_reads[pmax(1,f_inversions - num_inversions)]))==as.character(seqnames(inversions)), start(WW_reads[pmax(1,f_inversions - num_inversions)]),1)
        end_inversions <- ifelse(as.character(seqnames(WW_reads[pmin(length(WW_reads), p_inversions + num_inversions)]))==as.character(seqnames(inversions)),
                end(WW_reads[pmin(length(WW_reads), p_inversions + num_inversions)]), seqlengths(WW_reads)[as.character(seqnames(inversions))])

        #Assembling the new wide intervals
        wide_inversions <- GRanges(seqnames = seqnames(inversions), ranges=IRanges(start = start_inversions, end = end_inversions))


        #Finding reads that overlap the wide intervals
	#And then some annoying checks and formatting in case the wide intervals overlap 0 reads, because the pair2frgm step removes a few
        WC_select <- mergeByOverlaps(wide_inversions, WC_d)
	WC_missing <- wide_inversions[!overlapsAny(wide_inversions,WC_d)]
	WC_replace <- cbind(mergeByOverlaps(WC_missing, WC_missing, type="equal"), as.data.frame(matrix(20, ncol=6, nrow=length(WC_missing))))
        colnames(WC_replace) <- colnames(WC_select)
        WC_select <- rbind(WC_select, WC_replace)


        WW_select <- mergeByOverlaps(wide_inversions, WW_d)
	WW_missing <- inversions[!overlapsAny(wide_inversions,WW_d)]
	WW_replace <- cbind(mergeByOverlaps(WW_missing, WW_missing, type="equal"), as.data.frame(matrix(20, ncol=6, nrow=length(WW_missing))))
        colnames(WW_replace) <- colnames(WW_select)
        WW_select <- rbind(WW_select, WW_replace)

        #Finding the two highest deltaW peaks per interval
        WW_coords <- do.call(rbind, as.list(by(WW_select, WW_select[,1], function(x) boundaries(x))))
        WC_coords <- do.call(rbind, as.list(by(WC_select, WC_select[,1], function(x) boundaries(x))))

        if(genotype=="hom" | genotype=="low") {

                #The peaks give the new intervals
                new <-sort(GRanges(seqnames=WW_coords[,1], ranges=IRanges(start=WW_coords[,2], width=WW_coords[,3])))

        }

        if(genotype=="het") {

                WW_ranges <- sort(IRanges(start=WW_coords[,2], width=WW_coords[,3]))
                WC_ranges <- sort(IRanges(start=WC_coords[,2], width=WC_coords[,3]))


                intersection <- pintersect(WW_ranges, WC_ranges, resolve.empty="start.x")

                #If the intersection is short (<200 bp) then we just take the WW interval by default, since in general the WW intervals are larger.
                intersected_ranges <- IRanges(start=ifelse(width(intersection) < 200, start(WW_ranges), start(intersection)), width=ifelse(width(intersection) < 200,
                        width(WW_ranges),width(intersection)))

                new <-sort(GRanges(seqnames=WC_coords[,1], ranges=intersected_ranges))

        }

        #Removing new intervals if they don't overlap the original interval obtained by merging all overlapping inversions of the same genotype
        #A better way to do with would be to choose the two peaks so that the left one is left of end(inversions) and the right one is right of start(inversions)
        new <- new[pintersect(inversions,new)$hit]
        return(new)

}


############################################################################################################################################################################################

#Wrapper for the interval adjustment methods

adjust_intervals <- function(inversions, reads, confidence=0.95, base,  sex="female", prior, prior_male=c(0.5,0.5), adjust_method="deltas",paired_reads=TRUE, blacklist=""){ 

		#separating called inversions by genotype and probability and creating GRanges objects; I don't bother with reference homozygotes
		het <- inversions[inversions$probability>=confidence & inversions$probability!="blacklisted" & (inversions$genotype == "1|0" | inversions$genotype == "0|1"),]
		hom <- inversions[inversions$probability>=confidence & inversions$probability!="blacklisted" & (inversions$genotype == "1" | inversions$genotype == "1|1") ,]
		het <- GRanges(seqnames=het[,1], ranges=IRanges(start=het[,2], end=het[,3]), mcols= het[,-c(1:3)])
                hom <- GRanges(seqnames=hom[,1], ranges=IRanges(start=hom[,2], end=hom[,3]), mcols=hom[,-c(1:3)])

		if (adjust_method == "deltas" | adjust_method == "all" | adjust_method=="low") {
        	
			new_reads <- process_reads(reads[[1]],reads[[2]], paired_reads=paired_reads)
		}

		if (adjust_method == "deltas" | adjust_method == "all") {

		        new_het <- adjust_deltaW(new_reads[[1]], new_reads[[2]], new_reads[[3]], new_reads[[4]], inversions=het, genotype="het")
			new_hom <- adjust_deltaW(new_reads[[1]], new_reads[[2]], new_reads[[3]], new_reads[[4]], inversions=hom, genotype="hom")

			delta_inversions <- genotype_new_regions(list(new_het,new_hom), het, hom, reads, confidence, base, sex, prior, prior_male, blacklist=blacklist)		
			
			#Combining inversions whether or not the delta method worked. Either way it's worth trying to merge those that overlap.
			het <- delta_inversions[[1]]
			hom <- delta_inversions[[2]]				
	
		}

		#Need an empty dataframe for the last rbind, in case adjust_method isn't "all"
                new_low_ref <- data.frame(character(),integer(),integer(),integer(),integer(),integer(),integer(),character(),double())
                colnames(new_low_ref) <- c("chr","start","end",paste0(base[[2]],"_counts_W"), paste0(base[[2]],"_counts_C"), "WC_counts_W", "WC_counts_C","genotype","probability")


                if (adjust_method == "all" | adjust_method == "low" ) {

                        low <- inversions[inversions$probability<confidence & inversions$probability!="blacklisted",]
                        low <- GRanges(seqnames=low[,1], ranges=IRanges(start=low[,2], end=low[,3]), mcols= low[,-c(1:3)])
                        new_low <- adjust_deltaW(new_reads[[1]], new_reads[[2]], new_reads[[3]], new_reads[[4]], inversions=low, genotype="low")
                        new_low <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=new_low, background=base[[1]], base_state=base[[2]],  sex=sex, prior=prior, prior_male=prior_male)
                        new_low_het <- new_low[new_low$probability>=confidence & (new_low$genotype=="0|1" | new_low$genotype=="1|0"),]
                        new_low_hom <- new_low[new_low$probability>=confidence & (new_low$genotype=="1|1" | new_low$genotype=="1"),]
                        new_low_ref <- new_low[new_low$probability>=confidence & (new_low$genotype=="0|0" | new_low$genotype=="0"),]
		
	                new_low_het <- GRanges(seqnames=new_low_het[,1], ranges=IRanges(start=new_low_het[,2], end=new_low_het[,3]), mcols= new_low_het[,-c(1:3)])
			new_low_hom <- GRanges(seqnames=new_low_hom[,1], ranges=IRanges(start=new_low_hom[,2], end=new_low_hom[,3]), mcols= new_low_hom[,-c(1:3)])
			new_low_ref <- GRanges(seqnames=new_low_ref[,1], ranges=IRanges(start=new_low_ref[,2], end=new_low_ref[,3]), mcols= new_low_ref[,-c(1:3)])			
			
			het <- c(het,new_low_het)
			hom <- c(hom, new_low_hom)
                }

		if (adjust_method != "minimal") {

			#Merging inversions and checking whether they have the same genotype
			new_regions <- list(reduce(het), reduce(hom))
        	       	new_inversions <- genotype_new_regions(new_regions, het, hom, reads, confidence, base, sex, prior, prior_male, blacklist=blacklist)	
		
		} else {
	
			#Trying to find and genotype the intersection of all inversions in a cluster of overlaps
			new_regions <- list(minimal(het), minimal(hom))
                        new_inversions <- genotype_new_regions(new_regions, het, hom, reads, confidence, base, sex, prior, prior_male, blacklist=blacklist)
			
		}



		if (adjust_method == "deltas" | adjust_method == "all" | adjust_method == "low") {
			
			#If there are still overlapping inversions after trying adjustment and merging, just take the largest.
			het <- sort(new_inversions[[1]], by = ~width, decreasing=TRUE)
			hom <- sort(new_inversions[[2]], by = ~width, decreasing=TRUE)
			het_merged <- reduce(het)
			hom_merged <- reduce(hom)
			het <- het[findOverlaps(het_merged, het, select="first")]
                        hom <- hom[findOverlaps(hom_merged, hom, select="first")]
		
			#Combining the inversions from the different adjustment methods
			all <- c(het,hom)
		
		} else {

			#For "merge" and "minimal" don't bother about the ones that couldn't be merged
	                all <- c(new_inversions[[1]], new_inversions[[2]])

		}

		if (adjust_method == "low" | adjust_method == "all"){


			no_inversion <- inversions[inversions$probability<confidence | inversions$probability == "blacklisted" ,]
			no_inversion <- GRanges(seqnames=no_inversion[,1], ranges=IRanges(start=no_inversion[,2], end=no_inversion[,3]), mcols= no_inversion[,-c(1:3)])
			no_inversion <- no_inversion[is.na(findOverlaps(no_inversion, all, select="first")),]
			no_inversion <- no_inversion[is.na(findOverlaps(no_inversion, new_low_ref, select="first")),]
			no_inversion <- data.frame(chr=seqnames(no_inversion), start=start(no_inversion), end=end(no_inversion), mcols(no_inversion))
			new_low_ref <- data.frame(chr=seqnames(new_low_ref), start=start(new_low_ref), end=end(new_low_ref), mcols(new_low_ref))

		} else {
			
			no_inversion <- inversions[inversions$probability<confidence | inversions$probability == "blacklisted" ,]	

		}


		#Constructing a dataframe of the adjusted inversions and the reference homozygotes
               	all <- data.frame(chr=seqnames(all), start=start(all), end=end(all), mcols(all))
               	ref <- inversions[inversions$probability>=confidence & inversions$probability != "blacklisted" & ( inversions$genotype == "0|0" | inversions$genotype == "0"),]
		colnames(all) <- colnames(ref)
		colnames(no_inversion) <- colnames(ref)
		colnames(new_low_ref) <- colnames(ref)			
                inversions <- rbind(all, no_inversion, ref, new_low_ref)

                return(inversions)

        }

############################################################################################################################################################################################

#Checks whether adjusted inversions have the same genotype as the originals.

genotype_new_regions <- function(new_regions, het, hom, reads, confidence, base,  sex="female", prior, prior_male, paired_reads=TRUE, blacklist=""){

		#Re-genotyping the adjusted inversions
               	new_het <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=new_regions[[1]], background=base[[1]], base_state=base[[2]],  sex=sex,
                        prior=prior, prior_male=prior_male)
               	new_hom <- genotype_inversions(WW_reads=reads[[1]], WC_reads=reads[[2]], regions=new_regions[[2]], background=base[[1]], base_state=base[[2]],  sex=sex,
                        prior=prior, prior_male=prior_male)


		#Discarding the ones that don't match genotypes and constructing a dataframe
		new_het <- new_het[(new_het$genotype=="0|1" | new_het$genotype =="1|0") & new_het$probability >= confidence & new_het$probability != "blacklisted",]
		new_hom <- new_hom[(new_hom$genotype=="1|1" | new_hom$genotype =="1") & new_hom$probability >= confidence & new_hom$probability != "blacklisted",]
		new_het <- GRanges(seqnames=new_het[,1], ranges=IRanges(start=new_het[,2], end=new_het[,3]), mcols= new_het[,-c(1:3)])
		new_hom <- GRanges(seqnames=new_hom[,1], ranges=IRanges(start=new_hom[,2], end=new_hom[,3]), mcols=new_hom[,-c(1:3)])

		#Retrieving inversions that weren't successfully adjusted
		het <- het[is.na(findOverlaps(het, new_het, select="first")),]
                hom <- hom[is.na(findOverlaps(hom, new_hom, select="first")),]

		return(list(c(het,new_het), c(hom,new_hom)))

}

############################################################################################################################################################################################

#Given a vector with (possibly flat) peaks, this finds the two highest peaks
#It returns the coordinates of the right side of the leftmost peak and the left side of the rightmost peak
#Used to determine inversion start and end coordinates from a list of deltaw peaks
#Requires some spacing between the peaks (kind of messily)

twopeaks <- function(x){

t <- rle(x)

if(length(t$length) > 1){

	#retain only local maxima in t
	t$values[-(which(diff(sign(diff(t$values)))==-2)+1)] <- 0

	#Finding the highest peak
	m1 <- which.max(t$values)
	
	#Setting a rough spacer: peaks shouldn't be too close together!
	space <- as.integer(length(t$values)/12)
		
	#Set highest peak and nearby values to zero
	t$values[max(m1-space,1):min(m1+space,length(t$values))] <- 0

	#Finding the second-highest peak, and figuring out which one comes first
	m2 <- which.max(t$values)
	m <- sort(c(m1,m2))

	#Returning the inner coordinates of the peaks
	coords <- c(sum(t$lengths[1:m[1]]), sum(t$lengths[1:(m[2]-1)])+1)

} else if (length(t$length) == 1){

	coords <- (c(1, length(x))) 
}

return(coords)

}

############################################################################################################################################################################################

#Uses twopeaks to choose new start and end positions for an inversion

boundaries <- function(x){

#If there are no inversions, x$deltaW will have length zero; need to return an empty dataframe
if(length(x$deltaW)==0){

	return(data.frame(seqnames=character(), start=character(), width=character()))

}else if(length(x$deltaW)==1){

	return(data.frame(seqnames=as.character(seqnames(x[,1])),start=start(x[,1]),width=width(x[,1]) ))

}

peaks <- twopeaks(x$deltaW)

#Choosing input for IRanges based on twopeaks output
seqnames <- as.character(seqnames(x[,2][peaks[1]]))
start <- end(x[,2][peaks[1]])
end <- start(x[,2][peaks[2]])
width <- pmax(end - start + 1, 1)

return(data.frame(seqnames, start, width))

}


############################################################################################################################################################################################

#For a GRanges object containing multiple sets of overlapping intervals, finds the intersection of all intervals in each overlapping set, where possible

minimal <- function(interval){

#Group by overlaps
grouped <- findOverlaps(interval,reduce(interval),select="all")
grouped <- data.frame(queryHits(grouped), subjectHits(grouped))

#Separate sets of overlapping intervals and find the latest start and the earliest end
separated <- split(grouped,grouped[,2])
ranges <- cbind(sapply(c(1:length(separated)), function(x) as.character(seqnames(interval[separated[[x]][1,1]]))), sapply(c(1:length(separated)), function(x) max(start(interval[separated[[x]][,1]]))), sapply(c(1:length(separated)), function(x) min(end(interval[separated[[x]][,1]]))))

#If there is no minimal interval, replace it with an interval of length 0
ranges[as.numeric(ranges[,2]) >= as.numeric(ranges[,3]),3] <- 1
ranges[as.numeric(ranges[,2]) >= as.numeric(ranges[,3]),2] <- 1

#Creating a new set of intervals
new <- GRanges(seqnames=ranges[,1], ranges=IRanges(start=as.numeric(ranges[,2]),end=as.numeric(ranges[,3])))

return(new)

}


############################################################################################################################################################################################

#import a bed file

import <- function(path) {

bed <- read.table(path)
gr <- GRanges(seqnames=bed[,1], ranges=IRanges(start=bed[,2], end=bed[,3]))
return(gr)

}
