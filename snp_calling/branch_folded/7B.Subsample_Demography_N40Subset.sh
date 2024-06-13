#!/bin/bash

#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=10000mb
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00

# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/branch_folded/7B.Subsample_Demography_N40Subset.sh ${i}; done 
CHR=$1

#filtered vcfs directory
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs

outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs
mkdir $outvcfs $outvcfs/folded $outvcfs/unfolded $outvcfs/invariant_mm1 $outvcfs/raw $outvcfs/stats

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
bcftools view --threads 5 --samples-file ~/merondun/cuculus_migration/snp_calling/branch_folded/Samples_Demography_N10_CCW-CCE-COW-COE-Outgroups_2024MAR13.list -Ov ${vcfs}/${CHR}_raw.vcf.gz | \
        #RETAIN neutral sites 
        bedtools intersect -header -a - -b $neutral | \
        #EXCLUDE repeats
        bedtools subtract -header -a - -b $repeats | \
        #SET genotypes below MINDP to missing 
        bcftools +setGT -Ou -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #UPDATE AC fields after changing genotypes to missing 
        bcftools +fill-tags -Ou -- -t AC,AN  | \
        #RETAIN only sites with at least 90% genotypes
        bcftools view --threads 5 -i "MQ > 30 & F_MISSING < 0.1 & INFO/DP > ${low} & INFO/DP < ${high}" -Oz -o ${outvcfs}/raw/${CHR}_allsitesN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 ${outvcfs}/raw/${CHR}_allsitesN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### for FOLDED SFS, just subset 
bcftools view --threads 5 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -i 'QUAL > 20' -Oz -o ${outvcfs}/folded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/raw/${CHR}_allsitesN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 ${outvcfs}/folded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz


### for UNFOLDED SFS, polarize alleles based on C. poliocephalus
python ~/merondun/gen_gen/general/Polarize_VCF_Add_AA.py ~/merondun/cuculus_migration/snp_calling/Outgroups.list ${outvcfs}/folded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-A1tmp

# Switch reference and ancestral allele, ensure --recalc is added to switch AC fields
java -jar ~/modules/jvarkit/dist/vcffilterjdk.jar -f ~/modules/jvarkit/dist/script.js --recalc ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-A1tmp.vcf.gz --output ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-A2tmp.vcf.gz
bcftools index --threads 5 ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-A2tmp.vcf.gz

# Remove any alleles which aren't polarizable, and remove genotype fields which weren't switched (e.g. AD) because of ref allele switch
bcftools view --threads 5 --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 --max-af 0.999 -i 'INFO/AA!="U"' ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-A2tmp.vcf.gz -Ob -o - | \
    bcftools annotate -x ^FORMAT/GT,FORMAT/DP,FORMAT/GQ -o ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz
bcftools index --threads 5 ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz
#rm ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-*tmp.vcf.gz* 


### And then simply grab the invariant sites 
bcftools view --threads 5 --max-ac 0 -Oz -o $outvcfs/invariant_mm1/${CHR}_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/raw/${CHR}_allsitesN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 $outvcfs/invariant_mm1/${CHR}_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### Summarize SNP Counts in tidy format
raw=$(bcftools index -n ${vcfs}/${CHR}_raw.vcf.gz)
neutral=$(bcftools index -n ${outvcfs}/raw/${CHR}_allsitesN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)
folded=$(bcftools index -n ${outvcfs}/folded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)
unfolded=$(bcftools index -n ${outvcfs}/unfolded/${CHR}_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz)
invariant=$(bcftools index -n $outvcfs/invariant_mm1/${CHR}_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)

echo -e "${CHR}\tRaw\t${raw}" > ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tNeutral_MM1\t${neutral}" >> ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tFolded\t${folded}" >> ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tUnfolded\t${unfolded}" >> ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tInvariant\t${invariant}" >> ${outvcfs}/stats/${CHR}.snp.counts

