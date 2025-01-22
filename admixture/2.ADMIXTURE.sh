#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# submit as: for GROUP in $(cat GROUPS.list); do for K in {2..10}; do sbatch -J ADMIX_${GROUP}_${K} ~/merondun/cuculus_migration/admixture/2.ADMIXTURE.sh ${GROUP} ${K}; done ; done

GROUP=$1
K=$2

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/unrelated_n261/${GROUP}
bedfiles=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/autosomal_files

mkdir -p ${wd}
cd ${wd}

#Run Admixture
admixture -j7 --cv=5 ${bedfiles}/${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.bed ${K} > ${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ${bedfiles}/${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50 -fname ${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.${K}.P -qname ${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.${K}.Q -P 10 -o eval_${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50_${K}
