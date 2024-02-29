#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

#mamba activate momi-py36
SUFFIX=$1

pops=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/SampleSubset_DistanceK_n10_2023DEC06.pop
auto_vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/vcfs/Autosomes
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2/neutral__intergenic-4fold.bed
chr_map=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2/chr_map.txt

# Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 10 -Oz -o vcfs/Autosomes.${SUFFIX}.vcf.gz chromosome_vcfs/*.${SUFFIX}.vcf.gz
bcftools index --threads 10 vcfs/Autosomes.${SUFFIX}.vcf.gz 

#fix the chr_1 to simply '1', otherwise momi doesn't understand 
bcftools annotate --threads 10 --rename-chrs ${chr_map} -Oz -o vcfs/Autosomes.${SUFFIX}.tmp.vcf.gz vcfs/Autosomes.${SUFFIX}.vcf.gz
mv vcfs/Autosomes.${SUFFIX}.tmp.vcf.gz vcfs/Autosomes.${SUFFIX}.vcf.gz
bcftools index -f --threads 10 vcfs/Autosomes.${SUFFIX}.vcf.gz