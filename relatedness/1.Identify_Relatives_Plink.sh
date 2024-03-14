#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# mamba activate snps 
# to submit: for i in $(cat GROUPS.list); do sbatch -J REL_${i} ~/merondun/cuculus_migration/relatedness/1.Plink_Related.sh ${i} ; done 
GROUP=$1

raw_vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged
CHR=chr_10

awk '{OFS="\t"}{print "0", $1}' ${GROUP}.list > ${GROUP}.plink
#extract SNPs 
plink --keep ${GROUP}.plink --vcf $raw_vcfdir/${CHR}.vcf.gz --allow-extra-chr --const-fid --geno 0.1 --maf 0.05 --genome --min 0.185 --out ${GROUP}_relatedness

#and then merge the output files afterwards:
#mergem *miss | awk '{OFS="\t"}{print $1, $3, $5}' > ~/merondun/cuculus_migration/relatedness/Missingness_2024FEB29.txt
#mergem *rel | awk '{OFS="\t"}{print $1, $2, $7}' > ~/merondun/cuculus_migration/relatedness/Relatedness_2024FEB29.txt
