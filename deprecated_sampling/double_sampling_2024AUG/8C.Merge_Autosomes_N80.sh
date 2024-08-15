#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=3
#SBATCH --time=200:00:00

outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n80_vcfs
mkdir $outvcfs/merged

# FOLDED: Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 5 -Oz -o ${outvcfs}/merged/Autosomes_Folded_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/folded/*vcf.gz
bcftools index --threads 5 ${outvcfs}/merged/Autosomes_Folded_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

# INVARIANT: Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 5 -Oz -o ${outvcfs}/merged/Autosomes_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/invariant_mm1/*vcf.gz
bcftools index --threads 5 ${outvcfs}/merged/Autosomes_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz