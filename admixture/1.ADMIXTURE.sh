#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# for i in {1..10}; do sbatch -J ADMIXTURE_${i} ~/merondun/cuculus_migration/admixture/1.ADMIXTURE.sh ${i}; done
K=$1

cd admixture

#Run Admixture
admixture -j7 --cv=5 ../vcfs/Autosomes.MQ-5X-MM1-AA-LDr2w50.bed ${K} > Autosomes.MQ-5X-MM1-AA-LDr2w50.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../vcfs/Autosomes.MQ-5X-MM1-AA-LDr2w50 -fname Autosomes.MQ-5X-MM1-AA-LDr2w50.${K}.P -qname Autosomes.MQ-5X-MM1-AA-LDr2w50.${K}.Q -P 10 -o eval_Autosomes.MQ-5X-MM1-AA-LDr2w50_${K}