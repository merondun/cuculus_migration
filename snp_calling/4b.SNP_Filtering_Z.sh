#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/4.SNP_Filtering.sh ${i}; done 

CHR=$1

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses
cd ${WD}

mkdir -p chromosome_vcfs chromosome_vcfs/stats

raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged
chr_map=/dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/snp_calling/chr_map.txt

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#Subset samples
#bcftools view --threads 10 --samples-file ~/merondun/cuculus_migration/relatedness/Sample_List_Unrelated.list --force-samples -Ou ${raw_vcfs}/${CHR}.vcf.gz | \
#        #change chr_1 to 1 
#        bcftools annotate --threads 10 --rename-chrs ${chr_map} -Ou | \
#        bcftools norm -d both -Oz -o chromosome_vcfs/${CHR}_raw.vcf.gz 
#bcftools index --threads 10 chromosome_vcfs/${CHR}_raw.vcf.gz

#filter on DP, retain only sites within mean coverage + 2*sd or mean covearge - 2*sd 
#bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' chromosome_vcfs/${CHR}_raw.vcf.gz > chromosome_vcfs/${CHR}_raw.dp_stats.txt
#mean=$(cat chromosome_vcfs/${CHR}_raw.dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
#sd=$(cat chromosome_vcfs/${CHR}_raw.dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
#low=$(echo "$mean - 2*$sd" | bc)
#high=$(echo "$mean + 2*$sd" | bc)

#rm chromosome_vcfs/${CHR}_raw.dp_stats.txt

#subset SNPs, filter for MQ ONLY, fix chromosome name (stripping chr_), so new chromosome is simply 1 
#bcftools view --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 --threads 10 -Oz -o chromosome_vcfs/${CHR}_snp.MQ.vcf.gz -i "QUAL > 20 & MQ > 30 & INFO/DP > ${low} & INFO/DP < ${high}" chromosome_vcfs/${CHR}_raw.vcf.gz
#bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ.vcf.gz

#filter for >= 5x genotypes, apply 10% missingness filter
MINDP=3
bcftools +setGT -Ou chromosome_vcfs/${CHR}_snp.MQ.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #update AC fields after changing genotypes to missing 
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -Oz -o chromosome_vcfs/${CHR}_snp.MQ-3X-MM1.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-3X-MM1.vcf.gz

