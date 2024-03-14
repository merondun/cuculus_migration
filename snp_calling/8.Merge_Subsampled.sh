#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

#e.g. sbatch ~/merondun/cuculus_migration/snp_calling/8.Merge_Subsampled.sh MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased
SUFFIX=$1

#mamba activate momi-py36
vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs
autos=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/autosomal_files

# Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 20 -Oz -o $autos/Autosomes.${SUFFIX}.vcf.gz $vcfdir/*.${SUFFIX}.vcf.gz
bcftools index --threads 20 $autos/Autosomes.${SUFFIX}.vcf.gz 

# For admixture, do a MAF and LD prune 
bcftools view --threads 20 -Ou --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.999 $autos/Autosomes.${SUFFIX}.vcf.gz  | \
    #prune LD 
    bcftools +prune -m 0.2 --window 50 -Oz -o $autos/Autosomes.${SUFFIX}-LDr2w50.vcf.gz
bcftools index --threads 20 $autos/Autosomes.${SUFFIX}-LDr2w50.vcf.gz