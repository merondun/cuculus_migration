#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

CHR=chr_1
CHRnum=1
TMPDIR=/dss/dsslegfs01/pr53da/pr53da-dss-0021/tmp

vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs/unrelated_all/${CHR}_snp.MQ-5X.vcf.gz
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/species_distinction/no_pruning
list=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/species_distinction/ID.list

# No LD Pruning here 
bcftools view --threads 10 --samples-file ${list} -Ou ${vcf} | \
	bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -i 'F_MISSING < 0.2' -Oz -o ${wd}/${CHR}_snp.MQ-5X-NoPrune.vcf.gz
bcftools index --threads 10 ${wd}/${CHR}_snp.MQ-5X-NoPrune.vcf.gz

# Create plink bed files for admixture + perform pca 
plink --threads 10 --const-fid --vcf ${wd}/${CHR}_snp.MQ-5X-NoPrune.vcf.gz --chr-set 29 --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out ${wd}/${CHR}_snp.MQ-5X-NoPrune

# Calculate FST 
VCF=${wd}/chr_${CHRnum}_snp.MQ-5X-NoPrune.vcf.gz

length=$(bcftools view ${wd}/chr_${CHRnum}_snp.MQ-5X-NoPrune.vcf.gz | grep "ID=${CHRnum}," | sed 's/.*=//g' | sed 's/>//g')

tabix ${VCF}
pixy --stats fst --bypass_invariant_check yes --vcf ${VCF} --chunk_size ${length} --chromosomes ${CHRnum} --populations ID.pixypop --window_size ${length} --n_cores 10 --output_folder fst

cd ${wd}/admixture

for K in {2..10}; do 

echo "RUNNING K: ${K}"

#Run Admixture
admixture -j7 --cv=5 ../${CHR}_snp.MQ-5X-NoPrune.bed ${K} > ${CHR}_snp.MQ-5X-NoPrune.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../${CHR}_snp.MQ-5X-NoPrune -fname ${CHR}_snp.MQ-5X-NoPrune.${K}.P -qname ${CHR}_snp.MQ-5X-NoPrune.${K}.Q -P 10 -o eval_${CHR}_snp.MQ-5X-NoPrune.${K}

done 
