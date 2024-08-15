#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=3
#SBATCH --time=200:00:00
# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/branch_folded/7C.Subsample_Demography_N80Subset.sh ${i}; done 
CHR=$1

#filtered vcfs directory
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs

outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n80_vcfs
mkdir $outvcfs $outvcfs/folded $outvcfs/unfolded $outvcfs/invariant_mm1 $outvcfs/raw $outvcfs/stats

#samples 
samples=/dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/deprecated_sampling/double_sampling_2024AUG/Samples_Demography_N20_CCW-CCE-COW-COE_2024AUG12.list

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/GCA_017976375.1_bCucCan1.pri_genomic.CHR_strip.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

### Generate Depth statistics for filtering on DP, retain only sites within mean coverage + 2*sd or mean covearge - 2*sd 
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' ${vcfs}/${CHR}_raw.vcf.gz > ${outvcfs}/stats/${CHR}_raw.dp_stats.txt
mean=$(cat ${outvcfs}/stats/${CHR}_raw.dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
sd=$(cat ${outvcfs}/stats/${CHR}_raw.dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
# Round to nearest integer 
low=$(echo "$mean - 2*$sd" | bc) 
high=$(echo "$mean + 2*$sd" | bc)

#subset samples first, and filter for MQ. Right at the beginning, remove repeats and non-neutral sites 
MINDP=5
bcftools view --threads 5 --samples-file ${samples} -Ov ${vcfs}/${CHR}_raw.vcf.gz | \
        #RETAIN neutral sites 
        bedtools intersect -header -a - -b $neutral | \
        #EXCLUDE repeats
        bedtools subtract -header -a - -b $repeats | \
        #SET genotypes below MINDP to missing 
        bcftools +setGT -Ou -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #UPDATE AC fields after changing genotypes to missing 
        bcftools +fill-tags -Ou -- -t AC,AN  | \
        #RETAIN only sites with at least 90% genotypes
        bcftools view --threads 5 -i "MQ > 30 & F_MISSING < 0.1 & INFO/DP > ${low} & INFO/DP < ${high}" -Oz -o ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### for FOLDED SFS, just subset 
bcftools view --threads 5 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -i 'QUAL > 20' -Oz -o ${outvcfs}/folded/${CHR}_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 ${outvcfs}/folded/${CHR}_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### And then simply grab the invariant sites 
bcftools view --threads 5 --max-ac 0 -Oz -o $outvcfs/invariant_mm1/${CHR}_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 $outvcfs/invariant_mm1/${CHR}_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### Summarize SNP Counts in tidy format
raw=$(bcftools index -n ${vcfs}/${CHR}_raw.vcf.gz)
neutral=$(bcftools index -n ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)
folded=$(bcftools index -n ${outvcfs}/folded/${CHR}_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)
invariant=$(bcftools index -n $outvcfs/invariant_mm1/${CHR}_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)

echo -e "${CHR}\tRaw\t${raw}" > ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tNeutral_MM1\t${neutral}" >> ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tFolded\t${folded}" >> ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tInvariant\t${invariant}" >> ${outvcfs}/stats/${CHR}.snp.counts

