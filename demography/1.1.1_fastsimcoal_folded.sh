#!/bin/bash -l
#SBATCH -J WCcuckfolded
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/slurms/slurm-%j-%x.out

# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 B1 ${bs}; done 
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 W1 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 W2 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 W3 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 C1 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 C2 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 C3 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 B0 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 B2 ${bs}; done
# for bs in {1..50};do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 B3 ${bs}; done
# for i in {1..3}; do for bs in {21..50}; do sbatch 1.1.1_fastsimcoal_folded.sh CCW_CCE_COW_COE_folded 20,20,20,20 W${i}C ${bs}; done; done

conda activate easySFS

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics"

echo $(date)
STARTTIME=$(date +%s)

echo ${3}_run${bs}
cd ${dat}/sim_final
cd ${1}/fastsimcoal2

coalsim=1000000 #ideally between 200000 aand 1000000
ECM=100 # at least 20, better between 50 and 100

mkdir run${4}
cp ${3}.tpl ${3}.est ${3}_joint*AFpop*.obs run${4}
cd run${4}
${dat}/sim_final/fsc27 -t ${3}.tpl -e ${3}.est -m -M -L ${ECM} -n ${coalsim} -q -c 12 --foldedSFS


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
