#!/bin/bash -l
#SBATCH -J plot
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/slurms/slurm-%j-%x.out

# for i in {1..3}; do sbatch 1.2_plot.sh CCW_CCE_COW_COE_folded C${i} 50 CCW CCE COW COE; done
# for i in {0..1}; do sbatch 1.2_plot.sh CCW_CCE_COW_COE_folded B${i} 50 CCW CCE COW COE; done
# for i in {1..3}; do sbatch 1.2_plot.sh CCW_CCE_COW_COE_folded W${i}C 50 CCW CCE COW COE; done

conda activate easySFS

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/sim_final/${1}/fastsimcoal2

#create dummy files for runs which failed
for i in {1..50}
do
if [ $(ls ./run$i/${2}/$2.bestlhoods | wc -l) -eq 1 ]
then
    echo "file already exist"
else
    echo "generate dummy file"
    echo "9 0" > ./run$i/${2}/$2.bestlhoods
fi
done

#smallest difference between observed and expected likelihoods
cat run{1..100}/${2}/${2}.bestlhoods | grep -v "MaxObsLhood\|MaxEstLhood" | awk '{print NR,$(NF-1),$NF}' | awk '{print NR,$3-$2}' | sed -e "s/-9/NA/g" | grep -v "NA" | awk '{print $1, $2*1}' | sort -n -k2 > ${2}_${3}bs.bestlhoods
cat ./run{1..100}/${2}/${2}.bestlhoods | grep -v NPOP1 | nl | grep -v "9 0" > ${2}_${3}.bestlhoods.all

mkdir bestruns
for i in $(head -1 ${2}_${3}bs.bestlhoods | awk '{print $1}')
do
cp -r run${i}/${2} ./bestruns
cp run${i}/${2}.est ./bestruns/${2}
cp run${i}/${2}_*.obs run${i}/${2}.tpl run${i}/${2}.par ./bestruns/${2}
done

zless *.bestlhoods > ./bestruns/bestlhoods.txt

cd bestruns
conda activate r4.2
Rscript ${dat}/sim_final/calculateAIC.sh $2
Rscript ${dat}/sim_final/SFStools.R -t print2D -i $2 -z
Rscript ${dat}/sim_final/SFStools.R -t 2Dto1D -i $2 -z
Rscript ${dat}/sim_final/SFStools.R -t print1D -i $2 -z
Rscript ${dat}/sim_final/plotModel.R -p $2 -l ${4},${5},${6},${7}

zless ./*/*.AIC > AIC.txt
