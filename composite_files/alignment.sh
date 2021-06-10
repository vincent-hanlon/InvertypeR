# Removing adapters, trimming poor-quality bases, and excluding very short reads.
# TO DO: The user needs to make sure the adapter sequences are correct for any particular set of FASTQ files.

cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-o $1".trimmed.1.fastq.gz" \
        -p $1".trimmed.2.fastq.gz" \
	$1"_R1_001.fastq.gz" \
        $1"_R2_001.fastq.gz" \
	-m 30 \
	-q 15

######################################################################################################################
# Aligning reads to the GRCh38 reference genome.
# TO DO: The user needs to substitute the path to their copy of the bowtie-indexed reference genome.

bowtie2 -x '/projects/lansdorp/sspipe/refseq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index' -1 $1".trimmed.1.fastq.gz" -2 $1".trimmed.2.fastq.gz" -S $1".trimmed.sam"

#####################################################################################################################################
# Sorting SAM files and compressing them into BAM files (plus retaining only high-quality reads).

samtools view -q10 -f2 -bS $1".trimmed.sam" | samtools sort - -o $1".trimmed.sorted.bam"


#####################################################################################################################################
# Marking duplicate reads. 
# Need to include the path to your own jarfile.
# TO DO: The user needs to substitute the path to their copy of picard.jar.

mkdir ./metrics

java -jar '/home/vhanlon/programs/picard.jar' MarkDuplicates \
	I= $1".trimmed.sorted.bam" \
	O= $1".trimmed.sorted.mdup.bam" \
	M= metrics/$1".mdup.txt" \
	VALIDATION_STRINGENCY=LENIENT

######################################################################################################################################
# Indexing prepared BAM files.

samtools index $1".trimmed.sorted.mdup.bam"
