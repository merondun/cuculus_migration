#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

CHR=$1

WORKDIR=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop

cd ${WORKDIR}

mkdir full_vcf

raw_vcf_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/

#From the full gVCF, including invariant sites, subset only samples of interest
bcftools view --threads 5 --types snps --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${raw_vcf_dir}/${CHR}.SNPS.vcf.gz -Ou | \
       bcftools view --min-ac 1 --min-alleles 2 --max-alleles 2 -Oz -o full_vcf/${CHR}.raw.vcf.gz

#Apply filtering, only on MQ, DP
bcftools view --threads 5 --max-alleles 2 -i 'MQ > 30 && INFO/DP > 150 && QUAL > 20' -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz full_vcf/${CHR}.raw.vcf.gz
bcftools index full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz

#Phase VCF with beagle
java -jar -Xmx160g ~/modules/beagle.28Jun21.220.jar gt=full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz out=full_vcf/${CHR}.MQ20-DP150-Q20.PHASED nthreads=8 window=40 overlap=2 impute=true
bcftools index --threads 5 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz

#add the INFO DP,MQ and FMT/DP annotations back onto this VCF, from the pre-phased VCF
bcftools annotate --threads 5 -a full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz -c INFO/DP,INFO/MQ,FMT/DP full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
bcftools index --threads 5 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
