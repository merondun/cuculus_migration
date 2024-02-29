#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# mamba activate momi-py36
# submit as:   for i in $(cat chromosome_vcfs/SUFFIXES.list); do sbatch -J PREP_${i} ~/merondun/cuculus_migration/momi_demography/1.Prep_momi.sh ${i}; done
SUFFIX=$1

pops=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/SampleSubset_DistanceK_n10_2023DEC06.pop
auto_vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/vcfs/Autosomes
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2/neutral__intergenic-4fold.bed
chr_map=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2/chr_map.txt

#start with the merged autosomal VCF 
#prep momi file 
python -m momi.read_vcf --outgroup Outgroup --verbose ${auto_vcf}.${SUFFIX}.vcf.gz $pops ${outdir}/Autosomes.${SUFFIX}_counts.gz --bed ${neutral}

#create sfs file
python -m momi.extract_sfs ${outdir}/Autosomes.${SUFFIX}_sfs.gz 100 ${outdir}/Autosomes.${SUFFIX}_counts.gz

SNPS=$(bcftools index -n vcfs/Autosomes.${SUFFIX}.vcf.gz)
echo "FINISHED PREP FOR FILTER: ${SUFFIX} THERE ARE: ${SNPS} SNPS"