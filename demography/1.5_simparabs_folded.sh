#!/bin/bash -l
#SBATCH -J parabsfolded
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/slurms/slurm-%j-%x.out

# for i in {1..100}; do sbatch 1.5_simparabs_folded.sh W1x $i finite CCW_CCE_COW_COE_folded; done
# for i in {1..100}; do sbatch 1.5_simparabs_folded.sh B2 $i finite CCW_CCE_COW_COE_folded; done

conda activate easySFS

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/sim_final/${4}/fastsimcoal2
cd bestruns/modelfit
cd ${1}_${3}

cp ../${1}.tpl ./${1}_${3}_${2}/${1}_${3}.tpl
cp ../${1}.est ./${1}_${3}_${2}/${1}_${3}.est
cp ../${1}.pv ./${1}_${3}_${2}/${1}_${3}.pv

cd ${1}_${3}_${2}
echo ${1}_${3}_${2}
${dat}/sim_final/fsc27 -t ${1}_${3}.tpl -e ${1}_${3}.est –initvalues ${1}_${3}.pv -m -M -L 100 -n 1000000 -q -c 8 --foldedSFS 

##legend##
#-j: output one simulated or bootstrapped SFS per file in a separate directory for easier analysis
#-d: Computes the site frequency spectrum (SFS) for the derived alleles for each population sample and for all pairs of samples (joint 2D SFS)
#-s0: output all SNPs in the DNA sequence
#-x: Does not generate Arlequin output file
#-I:Generates DNA mutations according to an infinite site (IS) mutation model
#-q: quiet mode
#-n100: 100 simulations


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
