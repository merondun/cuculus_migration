#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/7.Subsample_Demography.sh ${i}; done 
CHR=$1

#filtered vcfs directory, polarized alleles
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs
outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs

mkdir $outvcfs

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

#Subset samples
bcftools view --threads 10 --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list -Ov ${vcfs}/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz | \
        #RETAIN neutral sites 
        bedtools intersect -header -a - -b $neutral | \
        #EXCLUDE repeats
        bedtools subtract -header -a - -b $repeats | \
        #re-apply 10% missing data threshold 
        bcftools view --threads 10 -i 'F_MISSING < 0.1' -Oz -o ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 10 ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz

#Phase VCF with beagle 
java -jar -Xmx160g ~/modules/beagle.28Jun21.220.jar gt=${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz out=${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp nthreads=8 window=40 overlap=2 impute=true
bcftools index --threads 10 ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz

#add the INFO DP,MQ and FMT/DP annotations back onto this VCF, from the pre-phased VCF 
bcftools annotate --threads 10 -a ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz -c INFO/DP,INFO/MQ,FMT/DP ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz -Oz -o ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.vcf.gz
bcftools index --threads 10 ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.vcf.gz

rm ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz.csi 