#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

#mamba activate snps 
GROUP=$1

raw_vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged
CHR=chr_1

bcftools view --threads 10 --force-samples --samples-file ${GROUP}.list $raw_vcfdir/${CHR}.vcf.gz -Ob | \
        bcftools view --types snps --threads 10 --min-af 0.3 --max-af 0.7 --min-alleles 2 --max-alleles 2 -i 'F_MISSING<0.1' -Ob | \
        bcftools +prune -m 0.1 -w 10kb --nsites-per-win 1 -o related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz
bcftools index --threads 10 related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz
bcftools index -n --threads 10 related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz

#extract phi and missingness rate 
vcftools --gzvcf related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz --relatedness2 --out related/${CHR}_${GROUP}
vcftools --gzvcf related/${CHR}_${GROUP}_MAF50-LDr1.vcf.gz --missing-indv --out related/${CHR}_${GROUP}

#format nicely for R 
awk -v g=${GROUP} '{OFS="\t"}{print $0, c, g}' related/${CHR}_${GROUP}.imiss > ${GROUP}.miss
awk -v g=${GROUP} '{OFS="\t"}{print $0, c, g}' related/${CHR}_${GROUP}.relatedness2 > ${GROUP}.rel

done 