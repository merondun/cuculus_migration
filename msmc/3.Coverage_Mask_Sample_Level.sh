#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

#submit sample as positional for RUN in $(cat ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list); do sbatch -J COV_${RUN} 2.Coverage_Masks.sh ${RUN} ; done
SAMPLE=$1

bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/Illumina_Alignments_Merged
#output mask directory
maskdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/coverage_masks
#directory with the n=40 subsampled phased, neutral and non-repetitive VCFS
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/full_vcf
#output directory with the individual chromosome level vcfs
indvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/individual_vcfs

cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop
mkdir $maskdir $maskdir/work $indvcfs

#loop through each chromosome
for CHR in $(cat Chromosomes.list); do

#create a coverage/MQ mask
mosdepth --threads 2 --chrom ${CHR} $maskdir/work/${SAMPLE}_${CHR} $bamdir/${SAMPLE}.bam
avg=$(awk '{print $4}' $maskdir/work/${SAMPLE}_${CHR}.mosdepth.summary.txt | tail -n 1)
min=$(printf "%.0f" $(echo "$avg / 2" | bc -l))
max=$(printf "%.0f" $(echo "$avg * 2" | bc -l))

echo "FOR SAMPLE: ${SAMPLE} average coverage is ${avg}, retainining positions between ${min} and ${max}"
#grab only sites greater than or below half the chromosome avg or double it. grab sites we want to RETAIN!
zcat $maskdir/work/${SAMPLE}_${CHR}.per-base.bed.gz | awk -v x=${max} -v n=${min} '$4 < x && $4 > n' | awk '{OFS="\t"}{print $1, $2, $3}' | bgzip -c > $maskdir/${SAMPLE}_${CHR}.bed.gz

#finally, subset that sample from the population VCF
bcftools view --threads 2 --samples ${SAMPLE} -Ob full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz | \
        bcftools view --threads 2 --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -Oz -o $indvcfs/${SAMPLE}_${CHR}.vcf.gz
bcftools index $indvcfs/${SAMPLE}_${CHR}.vcf.gz

done
