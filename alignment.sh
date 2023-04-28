# 1st argument ($1) is the number of CPUs for parallelizing
# 2nd argument ($2) is the path to the reference genome.

# Run like this:
# bash alignment.sh numCPU /path/to/reference/genome.fasta.

# At the moment, this is set up for PE reads, but it can easily be altered for SE reads.
# Similarly, this script expects FASTQ files ending in "_R1_001.fastq.gz" and "_R2_001.fastq.gz".

for i in *R1_001.fastq.gz
do
       name=$(echo $i | sed 's/_R1_001\.fastq\.gz//')
       bwa-mem2 mem -t "$1" "$2" "$i" "$name"_R2_001.fastq.gz > "$name".sam
done

# compresses, sorts, marks duplicates, and indexes files with aligned sequence reads.
process() {

        name=$(echo $1 | sed 's/\.sam//')

        samtools view -q10 -f2 -bS "$name".sam | samtools sort -n -o "$name".bam
        samtools fixmate -m "$name".bam - | samtools sort -o - - | samtools markdup - "$name".sorted.mdup.bam
        samtools index "$name".sorted.mdup.bam

}

export -f process

parallel --jobs "$1" process ::: *.sam
