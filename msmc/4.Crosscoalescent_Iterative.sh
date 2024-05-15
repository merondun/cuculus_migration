#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=120:00:00

#mamba activate samtools0.1.19
#submit with script 5, otherwise: sbatch ~/merondun/cuculus_migration/msmc/4.Crosscoalescent_Iterative.sh CCW CCE 1

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
#file with the mosdepth coverage masks (sites < 1/2 or > 2 the chromosome-sample-specific coverage are ignored )
maskdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/coverage_masks
#vcfs, individual for sample-chr
vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/individual_vcfs
#genome-wide mappability mask
gwmask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/demography/mappability/masks/
#neutral sites
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/neutral-sites_cuckoo__intergenic-intron-4fold.bed
#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin

#submit population 1, population 2, iteration
P1=$1
P2=$2
IT=$3

mkdir crosscoal crosscoal/input crosscoal/output

#grab 2 random samples from each population
grep -w ${P1} ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.pop | awk '{print $1}' | shuf | head -n 2 > crosscoal/${P1}_${P2}_${IT}.inds
grep -w ${P2} ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.pop | awk '{print $1}' | shuf | head -n 2 >> crosscoal/${P1}_${P2}_${IT}.inds

#loop through chromosomes and create the msmc input file
for i in $(cat Chromosomes.list); do

rm crosscoal/${P1}_${P2}_${IT}.${i}_POPmask.list crosscoal/${P1}_${P2}_${IT}.${i}_POPvcf.list

for j in $(cat crosscoal/${P1}_${P2}_${IT}.inds); do

echo "--mask=${maskdir}/${j}_${i}.bed.gz " >> crosscoal/${P1}_${P2}_${IT}.${i}_POPmask.list
echo "${vcfdir}/${j}_${i}.vcf.gz" >> crosscoal/${P1}_${P2}_${IT}.${i}_POPvcf.list

done

generate_multihetsep.py --chr ${i} --negative_mask $repeats --mask ${gwmask}/GCA_017976375.1_${i}.mask.bed.gz --mask $neutral $(cat crosscoal/${P1}_${P2}_${IT}.${i}_POPmask.list) $(cat crosscoal/${P1}_${P2}_${IT}.${i}_POPvcf.list) > crosscoal/input/${P1}_${P2}_${IT}_${i}.multihetsep.txt

done

#Run MSMC2, first on each population, then on both together (exhaustive haplotype comparisons)
~/modules/msmc2_Linux -t 15 -p 1*2+16*1+1*2 -I 0,1,2,3 -o crosscoal/output/${P1}_${P2}_${IT}_msmc_FIRST $(ls crosscoal/input/${P1}_${P2}_${IT}_*)
~/modules/msmc2_Linux -t 15 -p 1*2+16*1+1*2 -I 4,5,6,7 -o crosscoal/output/${P1}_${P2}_${IT}_msmc_SECOND $(ls crosscoal/input/${P1}_${P2}_${IT}_*)
~/modules/msmc2_Linux -t 15 -p 1*2+16*1+1*2 -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -s -o crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS $(ls crosscoal/input/${P1}_${P2}_${IT}_*)

#Merge all 4 iterations
python ~/modules/msmc-tools/combineCrossCoal.py \
    crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.final.txt \
    crosscoal/output/${P1}_${P2}_${IT}_msmc_FIRST.final.txt \
    crosscoal/output/${P1}_${P2}_${IT}_msmc_SECOND.final.txt > crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.CROSSCOAL_final.txt

#Also run MSMC-IM, sing default settings
mu=1.01e-08
python ~/modules/MSMC-IM/MSMC_IM.py crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.CROSSCOAL_final.txt -p 1*2+16*1+1*2 --printfittingdetails --plotfittingdetails --xlog -mu ${mu} -o crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.CROSSCOAL_MSMCIM.final.txt
