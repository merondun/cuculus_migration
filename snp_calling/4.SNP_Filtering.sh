#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/4.SNP_Filtering.sh ${i}; done 
CHR=$1

mkdir chromosome_vcfs chromosome_vcfs/stats

raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged
chr_map=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2/chr_map.txt

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#Subset samples
bcftools view --threads 10 --samples-file ~/merondun/cuculus_migration/relatedness/Sample_List_Unrelated.list --force-samples -Ou ${raw_vcfs}/${CHR}.vcf.gz | \
        #change chr_1 to 1 
        bcftools annotate --threads 10 --rename-chrs ${chr_map} -Ou | \
        bcftools norm -d both -Oz -o chromosome_vcfs/${CHR}_raw.vcf.gz 
bcftools index --threads 10 chromosome_vcfs/${CHR}_raw.vcf.gz

#filter on DP, retain only sites within mean coverage + 2*sd or mean covearge - 2*sd 
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' chromosome_vcfs/${CHR}_raw.vcf.gz > chromosome_vcfs/${CHR}_raw.dp_stats.txt
mean=$(cat chromosome_vcfs/${CHR}_raw.dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
sd=$(cat chromosome_vcfs/${CHR}_raw.dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
low=$(echo "$mean - 2*$sd" | bc)
high=$(echo "$mean + 2*$sd" | bc)

rm chromosome_vcfs/${CHR}_raw.dp_stats.txt

#subset SNPs, filter for MQ ONLY, fix chromosome name (stripping chr_), so new chromosome is simply 1 
bcftools view --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 --threads 10 -Oz -o chromosome_vcfs/${CHR}_snp.MQ.vcf.gz -i "QUAL > 20 & MQ > 30 & INFO/DP > ${low} & INFO/DP < ${high}" chromosome_vcfs/${CHR}_raw.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ.vcf.gz

#filter for >= 5x genotypes, apply 10% missingness filter
MINDP=5
bcftools +setGT -Ou chromosome_vcfs/${CHR}_snp.MQ.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #update AC fields after changing genotypes to missing 
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -i 'F_MISSING < 0.1' -Oz -o chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz

# Polarize alleles based on C. poliocephalus / C. micropterus shared derived alleles
python ~/merondun/gen_gen/general/Polarize_VCF_Add_AA.py ~/merondun/cuculus_migration/snp_calling/Outgroups.list chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A1

# Switch reference and ancestral allele, ensure --recalc is added to switch AC fields
java -jar ~/modules/jvarkit/dist/vcffilterjdk.jar -f ~/modules/jvarkit/dist/script.js --recalc chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A1.vcf.gz --output chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz

# Remove any alleles which aren't polarizable, and remove genotype fields which weren't switched (e.g. AD) because of ref allele switch
bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 --max-af 0.999 --threads 10 -i 'INFO/AA!="U"' chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz -Ob -o - | \
    bcftools annotate -x ^FORMAT/GT,FORMAT/DP,FORMAT/GQ -o chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz
rm chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A1.vcf.gz* chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz*

# LD Prune
bedtools intersect -header -a chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz -b $neutral | \
        bcftools +prune -m 0.2 --window 50 -Oz -o chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz
bcftools index --threads 10  chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz

# And phase with beagle 
java -Xmx40g -jar ~/modules/beagle.28Jun21.220.jar gt=chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz out=chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50-phased nthreads=10 impute=true
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50-phased.vcf.gz

# Grab invariant sites, filter for MQ, DP (same thresholds), and missingness
echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
bcftools view --max-ac 0 -i "MQ > 30 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" --threads 10 -Ou chromosome_vcfs/${CHR}_raw.vcf.gz -Oz -o chromosome_vcfs/${CHR}_1N.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_1N.vcf.gz

# re-merge the invariant and filtered SNPs
bcftools concat --threads 10 -Ob chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz chromosome_vcfs/${CHR}_1N.vcf.gz | \
        bcftools sort -Oz -o chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz

#SNP Counts in tidy format
raw=$(bcftools index -n chromosome_vcfs/${CHR}_raw.vcf.gz)
mq=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ.vcf.gz)
mm1=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz)
aa=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz)
ld=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz)
total=$(bcftools index -n chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz)

echo -e "${CHR}\tRaw\t${raw}" > chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tMQ20\t${mq}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\t5X_MM1\t${mm1}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tPolarized\t${aa}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tLD_R2_W50\t${ld}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tAllSites\t${total}" >> chromosome_vcfs/stats/${CHR}.snp.counts

#create simon's divergence input file
echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 2 --skipIndels -i chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz | \
        bgzip > chromosome_vcfs/${CHR}.geno.gz

rm chromosome_vcfs/${CHR}_1N.vcf.gz*