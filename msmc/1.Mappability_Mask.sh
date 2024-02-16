#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=100:00:00

progdir=/dss/dsshome1/lxc07/di39dux/modules/misc/seq/seqbility

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa


$progdir/splitfa ${genome} 50 | split -l 20000000
cat x* | gzip > simu_reads.fq.gz
bwa index ${genome}
bwa aln -R 1000000 -O 3 -E 3 ${genome} simu_reads.fq.gz > simu_reads.sai
bwa samse -f simu_reads.sam $genome simu_reads.sai simu_reads.fq.gz
perl $progdir/gen_raw_mask.pl simu_reads.sam > rawMask_50.fa
$progdir/gen_mask -l 50 -r 0.5 rawMask_50.fa > mask_50_50.fa
