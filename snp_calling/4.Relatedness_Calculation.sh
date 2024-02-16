#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

#chromosome (chr_10 here) 
i=$1

#loop by species, e.g. i=CC; therefore call relatives within CC, and CO separately 
for j in $(cat Groups.list); do 

#Subset tthat species, require MAF 5% 
bcftools view --threads 10 --force-samples --samples-file ${j}.list ../merged/${i}.vcf.gz -Ou | \
       bcftools view --threads 10 --min-af 0.05 --min-alleles 2 --max-alleles 2 -i 'F_MISSING<0.1' -Oz -o ${i}.${j}.vcf.gz
bcftools index ${i}.${j}.vcf.gz

#extract phi and missingness rate 
vcftools --gzvcf ${i}.${j}.vcf.gz --relatedness2 --out ${i}.${j}
vcftools --gzvcf ${i}.${j}.vcf.gz --missing-indv --out ${i}.${j}

#format nicely for R 
awk -v c=${i} -v g=${j} '{OFS="\t"}{print $0, c, g}' ${i}.${j}.imiss > ${i}.${j}.miss
awk -v c=${i} -v g=${j} '{OFS="\t"}{print $0, c, g}' ${i}.${j}.relatedness2 > ${i}.${j}.rel

done 