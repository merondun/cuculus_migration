#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

#submit sample as positional for RUN in $(cat ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list); do sbatch -J COV_${RUN} ~/merondun/cuculus_migration/msmc/3.Coverage_Mask_Sample_Level.sh ${RUN} ; done 
SAMPLE=$1

bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/Illumina_Alignments_Merged
#output mask directory
maskdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/coverage_masks
#directory with the n=40 subsampled phased, neutral and non-repetitive VCFS 
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs
#output directory with the individual chromosome level vcfs 
indvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/individual_vcfs

mkdir $maskdir $maskdir/work $indvcfs

#loop through each chromosome 
for CHR in $(cat Chromosomes.list); do 

#create a coverage/MQ mask
mosdepth --threads 5 --chrom ${CHR} $maskdir/work/${SAMPLE}_${CHR} $bamdir/${SAMPLE}.bam
avg=$(awk '{print $4}' $maskdir/work/${SAMPLE}_${CHR}.mosdepth.summary.txt | tail -n 1)
min=$(printf "%.0f" $(echo "$avg / 2" | bc -l))
max=$(printf "%.0f" $(echo "$avg * 2" | bc -l))

echo "FOR SAMPLE: ${SAMPLE} average coverage is ${avg}, retainining positions between ${min} and ${max}"
#grab only sites greater than or below half the chromosome avg or double it. grab sites we want to RETAIN! 
zcat $maskdir/work/${SAMPLE}_${CHR}.per-base.bed.gz | awk -v x=${max} -v n=${min} '$4 < x && $4 > n' | awk '{OFS="\t"}{print $1, $2, $3}' | bgzip -c > $maskdir/${SAMPLE}_${CHR}.bed.gz

#finally, subset that sample from the population VCF
bcftools view --threads 4 --samples ${SAMPLE} -Ob ${vcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.vcf.gz | \
        bcftools view --threads 4 --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -Oz -o $indvcfs/${SAMPLE}_${CHR}.vcf.gz
bcftools index $indvcfs/${SAMPLE}_${CHR}.vcf.gz

done 