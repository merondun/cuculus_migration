#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00
#SBATCH --output=/dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/pi_tajima/Slurms/slurm-%j.out

#mamba activate snps

# Directory with n=40subsample demography individuals. All sites, filtered, neutral, nonrepetitive, 10% max-missing
n40vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs

# Output directory 
divdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/divergence

# individuals $ID \t $population 
popfile=/dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.pop

mkdir $divdir $divdir/work $divdir/out

CHR=$1

#analyze in 25KB windows, requiring at least 10/25Kb present
sites=25000
missites=10000

echo "${CHR} working on ${sites} window with ${missites} missing sites max"


# Create geno file for popgen windows 
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing \
	        --ploidy 2 --skipIndels -i ${n40vcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz | \
		        bgzip > ${n40vcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.geno.gz

# Calculate pi/tajima 
popgenWindows.py -f phased -g ${n40vcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.geno.gz \
	        -o $divdir/work/${CHR}.csv.gz --analysis popFreq popDist popPairDist \
		        --windType coordinate -w $sites -m $missites --ploidy 2 -T 5 -p CCW -p CCE -p COW -p COE --popsFile ${popfile}

# Just convert to tab
zcat $divdir/work/${CHR}.csv.gz | tr ',' '\t' > $divdir/out/${CHR}.txt
