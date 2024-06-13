#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

# mamba activate snps

#filtered vcfs directory
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs

outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs
mkdir $outvcfs $outvcfs/folded $outvcfs/unfolded $outvcfs/invariant_mm1 $outvcfs/raw $outvcfs/stats

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/GCA_017976375.1_bCucCan1.pri_genomic.CHR_strip.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

### Subsample N40, folded
bcftools view --threads 5 --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${outvcfs}/merged/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz | \
    bcftools view --threads 5 --min-ac 1 --max-af 0.9999 -Oz -o $outvcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 $outvcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### Grab the invariant sites within the N=40 
bcftools view --threads 5 --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${outvcfs}/merged/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz | \
    bcftools view --threads 5 --max-ac 0 -Oz -o $outvcfs/merged/n40/Autosomes_invariantN40-SNPOutgroup.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz 
bcftools index --threads 5 $outvcfs/merged/n40/Autosomes_invariantN40-SNPOutgroup.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### Subsample N40, unfolded 
bcftools view --threads 5 --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${outvcfs}/merged/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz | \
    bcftools view --threads 5 --min-ac 1 --max-af 0.9999 -Oz -o $outvcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz
bcftools index --threads 5 $outvcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz

### For admixture, do a MAF and LD prune
bcftools view --threads 5 -Ou --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 $outvcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz | \
        #prune LD
    bcftools +prune -m 0.2 --window 50 -Oz -o ${outvcfs}/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-MAF5-LDr2w50.vcf.gz
bcftools index --threads 5 ${outvcfs}/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-MAF5-LDr2w50.vcf.gz

### For admixture, do a MAF and LD prune, unfolded
bcftools view --threads 5 -Ou --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.95 $outvcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz | \
        #prune LD
    bcftools +prune -m 0.2 --window 50 -Oz -o $outvcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA-MAF5-LDr2w50.vcf.gz
bcftools index --threads 5 $outvcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA-MAF5-LDr2w50.vcf.gz

### Create plink v1 bed files for admixture + perform pca 
plink --threads 5 --vcf ${outvcfs}/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-MAF5-LDr2w50.vcf.gz --chr-set 29 --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out ${outvcfs}/merged/n40/plink/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-MAF5-LDr2w50

### Create plink v1 bed files for admixture + perform pca, unfolded
plink --threads 5 --vcf $outvcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA-MAF5-LDr2w50.vcf.gz --chr-set 29 --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out $outvcfs/merged/n40/plink/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA-MAF5-LDr2w50

### Grab the invariant sites within the N=40 
bcftools view --threads 5 --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${outvcfs}/merged/Autosomes_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz | \
    bcftools view --threads 5 --max-ac 0 -Oz -o $outvcfs/merged/n40/Autosomes_invariantN40-WithOG.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 $outvcfs/merged/n40/Autosomes_invariantN40-WithOG.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### Merge invariant into one file
bcftools concat --threads 5 -a -Ou $outvcfs/merged/n40/Autosomes_invariantN40-SNPOutgroup.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz $outvcfs/merged/n40/Autosomes_invariantN40-WithOG.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz | \
        bcftools sort -Oz -o $outvcfs/merged/n40/Autosomes_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 $outvcfs/merged/n40/Autosomes_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz