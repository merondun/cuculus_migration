#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# submit as:  for GROUP in $(cat GROUPS.list); do sbatch -J MERGE_${i} ~/merondun/cuculus_migration/admixture/1.Subset_Subgroups.sh ${GROUP}; done 
# mamba activate snps

GROUP=$1

mkdir autosomal_files 

#intergenic and 4-fold degenerate sites 
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

# Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 20 -Ou chromosome_vcfs/*AA.vcf.gz | \
    #subset samples 
    bcftools view --threads 20 -Ou --force-samples --samples-file ~/merondun/cuculus_migration/relatedness/${GROUP}_Unrelated.list | \
    #ensure only variant SNPs are left 
    bcftools view --threads 20 -Ov --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.999 -i 'F_MISSING < 0.1' | \
    #intersect with neutral sites 
    bedtools intersect -header -a - -b $neutral | \
    #prune LD 
    bcftools +prune -m 0.2 --window 50 -Oz -o autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.vcf.gz
bcftools index --threads 20 autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.vcf.gz

#Create plink bed files for admixture + perform pca 
~/modules/plink2 --threads 20 --vcf autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.vcf.gz --chr-set 29 --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50