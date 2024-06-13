#!/bin/bash

#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=15000mb
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00

outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs

# FOLDED: Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 5 -Oz -o ${outvcfs}/merged/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/folded/*vcf.gz
bcftools index --threads 5  ${outvcfs}/merged/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz


# UNFOLDED: Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 5 -Oz -o ${outvcfs}/merged/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz ${outvcfs}/unfolded/*vcf.gz
bcftools index --threads 5  ${outvcfs}/merged/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz


# INVARIANT: Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 5 -Oz -o ${outvcfs}/merged/Autosomes_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/invariant_mm1/*vcf.gz
bcftools index --threads 5 ${outvcfs}/merged/Autosomes_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

# For admixture, do a MAF and LD prune
bcftools view --threads 5 -Ou --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 ${outvcfs}/merged/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz | \
	    #prune LD
    bcftools +prune -m 0.2 --window 50 -Oz -o ${outvcfs}/merged/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-LDr2w50.vcf.gz

