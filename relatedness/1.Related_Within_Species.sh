#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# mamba activate snps 
# to submit: for i in $(cat GROUPS.list); do sbatch -J REL_${i} ~/merondun/cuculus_migration/relatedness/1.Related_Within_Species.sh ${i} ; done 
GROUP=$1

raw_vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged
CHR=chr_10

mkdir related vcfs 

MINDP=10
#extract SNPs 
bcftools view --threads 20 --force-samples --samples-file ${GROUP}.list $raw_vcfdir/${CHR}.vcf.gz -Ou | \
        bcftools view --threads 20 --types snps --min-ac 1 --min-alleles 2 --max-alleles 2 -Oz -o vcfs/${GROUP}.vcf.gz
bcftools +setGT -Ou vcfs/${GROUP}.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --threads 20 --types snps --min-af 0.01 --min-alleles 2 --max-alleles 2 -i 'F_MISSING<0.1' -Oz -o vcfs/${GROUP}_MM1.vcf.gz
bcftools index --threads 5 vcfs/${GROUP}_MM1.vcf.gz
bcftools index -n --threads 5 vcfs/${GROUP}_MM1.vcf.gz

#extract phi and missingness rate 
vcftools --gzvcf vcfs/${GROUP}_MM1.vcf.gz --relatedness2 --out related/${GROUP}
vcftools --gzvcf vcfs/${GROUP}_MM1.vcf.gz --missing-indv --out related/${GROUP}

#format nicely for R 
awk -v g=${GROUP} '{OFS="\t"}{print $0, g}' related/${GROUP}.imiss > related/${GROUP}.miss
awk -v g=${GROUP} '{OFS="\t"}{print $0, g}' related/${GROUP}.relatedness2 > related/${GROUP}.rel


####afterwards, merge them:
# mergem *rel | awk '{OFS="\t"}{print $1, $2, $7}' > ~/merondun/cuculus_migration/relatedness/Relatedness_2024FEB29.txt
# mergem *miss | awk '{OFS="\t"}{print $1, $3, $5}' > ~/merondun/cuculus_migration/relatedness/Missingness_2024FEB29.txt