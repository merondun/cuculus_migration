#!/bin/bash -l
#SBATCH -J modelfit
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/slurms/slurm-%j-%x.out

# sbatch 1.4_modelfit_folded.sh CCW_CCE_COW_COE_folded W1
# sbatch 1.4_modelfit_folded.sh CCW_CCE_COW_COE_folded W1x
# sbatch 1.4_modelfit_folded.sh CCW_CCE_COW_COE_folded B2

conda activate easySFS

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/sim_final/${1}/fastsimcoal2
cd bestruns/modelfit

#infinite sites
${dat}/sim_final/fsc27 -i ${2}_infinite.par -j -m -s0 -x -I -q -n100 -c 2

#finite sites
${dat}/sim_final/fsc27 -i ${2}_finite.par -j -m -s0 -x -q -n100 -c 2

##legend##
#-j: output one simulated or bootstrapped SFS per file in a separate directory for easier analysis
#-m: Computes the site frequency spectrum (SFS) for th.pare minor alleles for each population sample and for all pairs of samples (joint SFS)
#-s0: output all SNPs in the DNA sequence
#-x: Does not generate Arlequin output file
#-I:Generates DNA mutations according to an infinite site (IS) mutation model
#-q: quiet mode
#-n100: 100 simulations

#remove the "no. of removed sites" text on finite SFS
cd $2_finite 
for i in {1..100}
do
cd $2_finite_${i}
sed -i 's/(.*)//g' $2_finite_jointMAFpop1_0.obs
sed -i 's/(.*)//g' $2_finite_jointMAFpop2_0.obs
sed -i 's/(.*)//g' $2_finite_jointMAFpop3_0.obs
sed -i 's/(.*)//g' $2_finite_jointMAFpop2_1.obs
sed -i 's/(.*)//g' $2_finite_jointMAFpop3_1.obs
sed -i 's/(.*)//g' $2_finite_jointMAFpop3_2.obs
cd ../
done

cd ${dat}/sim_final/${1}/fastsimcoal2/bestruns/modelfit
conda activate r4.2

for site in finite infinite
do 
cd $2_${site}
awk 'NR>1' $2_${site}.lhoodObs | tr '\t' '\n' | awk 'NR>2' | nl | sort -n -k2 -r > $2_${site}.lhoodObs_sorted
for i in $(head -1 $2_${site}.lhoodObs_sorted | awk '{print $1}')
do
mkdir ${2}
cp *_${i}/*obs ./${2}
cd ./${2}
rename ${2}_${site}_ ${2}_ *.obs
rename .obs .txt *.obs
cd ../
cp ../${2}*.obs ./${2}
Rscript ${dat}/sim_final/SFStools.R -t print1D -i $2 -z
done
cd ../
done
