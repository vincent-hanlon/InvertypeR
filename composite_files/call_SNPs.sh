ls ../*.bam > samples.list

callvariants.sh list=samples.list  ref=$1 out=called_snps.vcf ploidy=2 calldel=f callins=f
