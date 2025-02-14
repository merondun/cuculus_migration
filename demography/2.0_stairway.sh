#!/bin/bash -l
#SBATCH -J stairway
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/slurms/slurm-%j-%x.out

# sbatch 2.0_stairway.sh CCE
# sbatch 2.0_stairway.sh CCW
# sbatch 2.0_stairway.sh COE
# sbatch 2.0_stairway.sh COW
# sbatch 2.0_stairway.sh CO
# sbatch 2.0_stairway.sh CC
# sbatch 2.0_stairway.sh CCE_n20
# sbatch 2.0_stairway.sh CCW_n20
# sbatch 2.0_stairway.sh COE_n20
# sbatch 2.0_stairway.sh COW_n20

conda activate biotools 
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics"

cd $dat/sim_final/stairway

java -cp stairway_plot_es Stairbuilder ${1}_folded.blueprint
bash ${1}_folded.blueprint.sh
bash ${1}_folded.blueprint.plot.sh

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"

