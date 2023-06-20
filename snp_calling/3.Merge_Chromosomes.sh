#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir ../merged
CHR=$1
echo "${CHR}"
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/${CHR}.consmask
bcftools concat --file-list VCFs.list --threads 5 -a -r ${CHR} -Ou | \
        bcftools sort -Oz -o ../merged/${CHR}.vcf.gz -
bcftools index --threads 5 ../merged/${CHR}.vcf.gz
