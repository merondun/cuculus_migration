#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# Merge the VCFs into a single autosomal SNP set
bcftools concat --file-list Autosomes.list --threads 20 -Oz -o vcfs/Autosomes.IF-GF-MM1-AA-BP.vcf.gz
bcftools index --threads 20 vcfs/Autosomes.IF-GF-MM1-AA-BP.vcf.gz

VCF=vcfs/Autosomes.IF-GF-MM1-AA-BP.vcf.gz
#neutral sites 
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/neutral/neutral_cuckoo__intergenic-intron-4fold.bed

bedtools intersect -header -a ${VCF} -b $neutral | bcftools norm -d both -Oz -o vcfs/Autosomes.IF-GF-MM1-AA-BP-NEUTRAL.vcf.gz
bcftools index vcfs/Autosomes.IF-GF-MM1-AA-BP-NEUTRAL.vcf.gz
