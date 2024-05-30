#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=/dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/Slurms/slurm-%j.out

# mamba activate snps
# Submit as: for CHR in $(cat ~/merondun/cuculus_migration/Chromosomes.list); do sbatch -J SUB_${CHR} 3B.Create_Invariant_Bed_Subset_Samples.sh ${CHR}; done

CHR=$1

subvcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs
msmcdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin_populationlevel

echo "Creating invariant sites for: ${CHR}"

# Grab invariant sites 
bcftools view --threads 2 --max-ac 1 -Oz ${subvcf}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.vcf.gz | \
	bcftools query -f '%CHROM\t%POS\t%POS\n' | \
	# Convert VCF 1-based to bed 0-based
	awk '{OFS="\t"}{print $1, $2-1, $2}' | \
	# Merge overlapping intervals 
	bedtools merge -i - | bgzip -c > ${msmcdir}/vcfs/${CHR}.invariant.bed.gz
tabix ${msmcdir}/vcfs/${CHR}.invariant.bed.gz

for SAMPLE in $(cat ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list); do 

	echo "Subsetting VCF for sample: ${SAMPLE}"

	# Subset each sample into a SNP file 
	bcftools view --threads 2 --samples ${SAMPLE} -Ob ${subvcf}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.vcf.gz | \
	        bcftools view --threads 2 --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -Oz -o ${msmcdir}/vcfs/${SAMPLE}_${CHR}.vcf.gz
	bcftools index ${msmcdir}/vcfs/${SAMPLE}_${CHR}.vcf.gz

done 
