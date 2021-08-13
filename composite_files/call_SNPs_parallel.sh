#Calling SNPs with minimal parameters for shallow sequence data


ls ../*bam > samples.list

freebayes-parallel <(fasta_generate_regions.py $1.fai 5000000 |  grep -v 'chrUn\|random\|EBV\|chrM') $2 -f $1 -n 2 --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous -L samples.list | \
	bcftools view -H -v snps -m2 -M2 -i 'FORMAT/AD[0:0] > 0 && FORMAT/AD[0:1] > 1' > snps.vcf
