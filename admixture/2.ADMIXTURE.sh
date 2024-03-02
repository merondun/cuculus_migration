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

mkdir admixture admixture/${GROUP}
cd admixture/${GROUP}

#Run Admixture
admixture -j7 --cv=5 ../../autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.bed ${K} > ${GROUP}.MQ-5X-MM1-AA-LDr2w50.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../../autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50 -fname ${GROUP}.MQ-5X-MM1-AA-LDr2w50.${K}.P -qname ${GROUP}.MQ-5X-MM1-AA-LDr2w50.${K}.Q -P 10 -o eval_${GROUP}.MQ-5X-MM1-AA-LDr2w50_${K}