!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

CHR=$1

raw_vcf_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged

#From the full gVCF, including invariant sites, subset only samples of interest
bcftools view --threads 5 --samples-file Samples.list ${raw_vcf_dir}/chr_19.vcf.gz

#Apply filtering, only on MQ, DP
bcftools view --threads 5 --max-alleles 2 -i 'MQ > 20 && INFO/DP > 150 && QUAL > 20' -Oz -o full_vcf/chr_19.MQ20-DP150-Q20.vcf.gz

#Phase VCF with beagle 
java -jar -Xmx160g ~/modules/beagle.28Jun21.220.jar gt=full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz out=full_vcf/${CHR}.MQ20-DP150-Q20.PHASED nthreads=8 window=40 overlap=2 impute=true
bcftools index --threads 10 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz

#add the INFO DP,MQ and FMT/DP annotations back onto this VCF, from the pre-phased VCF 
bcftools annotate --threads 10 -a full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz -c INFO/DP,INFO/MQ,FMT/DP full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
bcftools index --threads 10 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz