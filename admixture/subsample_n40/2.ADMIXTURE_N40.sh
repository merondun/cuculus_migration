#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# mamba activate snps
# submit as: for RUN in $(cat RUNS.list); do for K in {2..10}; do sbatch -J ADMIX_${RUN}_${K} ~/merondun/cuculus_migration/admixture/subsample_n40/2.ADMIXTURE_N40.sh ${RUN} ${K}; done; done

RUN=$1
K=$2

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/unrelated_chyiyin_n40
plinkdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs/merged/n40/plink

cd $wd

mkdir ${RUN}
cd ${RUN}

#Run Admixture
admixture -j7 --cv=5 ${plinkdir}/${RUN}.bed ${K} > ${RUN}.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ${plinkdir}/${RUN} -fname ${RUN}.${K}.P -qname ${RUN}.${K}.Q -P 10 -o eval_${RUN}_${K}