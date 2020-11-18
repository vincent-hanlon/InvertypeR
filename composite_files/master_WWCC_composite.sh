#Master script to create a WW or CC composite BAM file from a set of strand-seq libraries for a sample

################################## User sets variables here #############################################

#Set the variables below as needed

threads=12
paired="TRUE"
scripts="/projects/lansdorp/composite_BAM_scripts/"
ref="/projects/lansdorp/sspipe/refseq/GRCh38.fasta"
sex="male"

##########################################################################################################
##########################################################################################################

Rscript $scripts/bpr.R 20000000 50 $threads $paired $sex "WW" $scripts

Rscript $scripts/printWWCCregions.R

cat <(awk '{print $4,$1,$2,$3,"WW"}' OFS="\t" ww_regions.txt) <(awk '{print $4,$1,$2,$3,"CC"}' OFS="\t" cc_regions.txt) > wwcc_regions.txt

#Merging WW and CC regions of BAM files, then reversing all reads in the CC BAM file to merge it with the WW reads
bash $scripts/WWCC_merge.sh wwcc_regions.txt $threads $paired


mkdir WW_CC
mv merged.WW.CC.bam* WW_CC

cd WW_CC

Rscript $scripts/bpr.R 100000 50 $threads $paired $sex "WW" $scripts
