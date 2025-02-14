#!/bin/bash -l
#SBATCH -J fsimplot
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/slurms/slurm-%j-%x.out

# sbatch 1.6_simparabs_plot.sh CCW_CCE_COW_COE_folded W1x finite

conda activate easySFS

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics"

echo $(date)
STARTTIME=$(date +%s)

popmodel=${2}
popgrp=${popmodel:0:4}
echo ${popgrp}

cd ${dat}/sim_final/${1}/fastsimcoal2/bestruns/modelfit/${2}_${3}

#create dummy files for runs which failed
for i in {1..100}
do
if [ $(ls ./${2}_${3}_${i}/${2}_${3}/${2}_${3}.bestlhoods | wc -l) -eq 1 ]
then
    echo "file already exist"
else
    echo "generate dummy file"
    mkdir ./${2}_${3}_${i}/${2}_${3}
    echo "9 0" > ./${2}_${3}_${i}/${2}_${3}/${2}_${3}.bestlhoods
fi
done

#smallest difference between observed and expected likelihoods
cat ${2}_${3}_{1..100}/${2}_${3}/${2}_${3}.bestlhoods | grep -v "MaxObsLhood\|MaxEstLhood" | awk '{print NR,$(NF-1),$NF}' | awk '{print NR,$3-$2}' | sed -e "s/-9/NA/g" | grep -v "NA" | sort -n -k2 > ${2}_${3}.bestlhoods
cat ./${2}_${3}_{1..100}/${2}_${3}/${2}_${3}.bestlhoods | grep -v NPOP1 | nl | grep -v "9 0" > ${2}_${3}.bestlhoods.all

