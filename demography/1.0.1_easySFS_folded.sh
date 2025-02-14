#!/bin/bash -l
#SBATCH -J foldedSFS
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/slurms/slurm-%j-%x.out

# sbatch 1.0.1_easySFS_folded.sh analyses/subsampled_n40_vcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive CCW_CCE_COW_COE_folded2 20,20,20,20 B1
# sbatch 1.0.1_easySFS_folded.sh analyses/subsampled_n40_vcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive CCW_CCE_CO_folded 20,20,40 B0
# sbatch 1.0.1_easySFS_folded.sh analyses/subsampled_n40_vcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive CC_CO_folded 40,40 B0
# sbatch 1.0.1_easySFS_folded.sh analyses/subsampled_n80_vcfs/merged/Autosomes_Folded_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive CCW_CCE_COW_COE_n80_folded 40,40,40,40 B0

conda activate easySFS

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/sim_final

#if folded SFS
python3 easySFS.py -i ${dat}/${1}.vcf.gz -p ${2}.pop -a -f --proj ${3} -o $2 --prefix ${4}

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
