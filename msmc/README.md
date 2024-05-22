# Coalescent analyses

1. Prep genome, create mappability mask
1A. Python script required for mappability mask
2. Using the all-sites raw VCF, subset samples, filter, and phase
3. Create a sample and chromosome-specific coverage mask file (with retained sites)
4. Run crosscoalescent and MSMC-IM, iteratively
5. Submit script 4 using all the population pairs across iterations

### Subset Samples

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

CHR=$1

mkdir full_vcf

raw_vcf_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/

#From the full gVCF, including invariant sites, subset only samples of interest
bcftools view --threads 10 --types snps --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${raw_vcf_dir}/${CHR}.SNPS.vcf.gz -Ou | \
       bcftools view --min-ac 1 --min-alleles 2 --max-alleles 2 -Oz -o full_vcf/${CHR}.raw.vcf.gz

#Apply filtering, only on MQ, DP
bcftools view --threads 10 --max-alleles 2 -i 'MQ > 30 && INFO/DP > 150 && QUAL > 20' -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz full_vcf/${CHR}.raw.vcf.gz
bcftools index full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz

#Phase VCF with beagle
java -jar -Xmx160g ~/modules/beagle.28Jun21.220.jar gt=full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz out=full_vcf/${CHR}.MQ20-DP150-Q20.PHASED nthreads=8 window=40 overlap=2 impute=true
bcftools index --threads 10 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz

#add the INFO DP,MQ and FMT/DP annotations back onto this VCF, from the pre-phased VCF
bcftools annotate --threads 10 -a full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz -c INFO/DP,INFO/MQ,FMT/DP full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
bcftools index --threads 10 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
```

### Create MSMC2 Masks

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

#submit sample as positional for RUN in $(cat ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list); do sbatch -J COV_${RUN} 2.Coverage_Masks.sh ${RUN} ; done
SAMPLE=$1

bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/Illumina_Alignments_Merged
#output mask directory
maskdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/coverage_masks
#directory with the n=40 subsampled phased, neutral and non-repetitive VCFS
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/full_vcf
#output directory with the individual chromosome level vcfs
indvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/individual_vcfs

cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin
mkdir $maskdir $maskdir/work $indvcfs

#loop through each chromosome
for CHR in $(cat Chromosomes.list); do

#create a coverage/MQ mask
mosdepth --threads 2 --chrom ${CHR} $maskdir/work/${SAMPLE}_${CHR} $bamdir/${SAMPLE}.bam
avg=$(awk '{print $4}' $maskdir/work/${SAMPLE}_${CHR}.mosdepth.summary.txt | tail -n 1)
min=$(printf "%.0f" $(echo "$avg / 2" | bc -l))
max=$(printf "%.0f" $(echo "$avg * 2" | bc -l))

echo "FOR SAMPLE: ${SAMPLE} average coverage is ${avg}, retainining positions between ${min} and ${max}"
#grab only sites greater than or below half the chromosome avg or double it. grab sites we want to RETAIN!
zcat $maskdir/work/${SAMPLE}_${CHR}.per-base.bed.gz | awk -v x=${max} -v n=${min} '$4 < x && $4 > n' | awk '{OFS="\t"}{print $1, $2, $3}' | bgzip -c > $maskdir/${SAMPLE}_${CHR}.bed.gz

#finally, subset that sample from the population VCF
bcftools view --threads 2 --samples ${SAMPLE} -Ob full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz | \
        bcftools view --threads 2 --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -Oz -o $indvcfs/${SAMPLE}_${CHR}.vcf.gz
bcftools index $indvcfs/${SAMPLE}_${CHR}.vcf.gz

done
```

### Run MSMC2

```bash
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

```

Submit with this:

```bash
# Loop through each pair of populations
for pair in "CCW CCE" "CCW COW" "CCW COE" "CCE COW" "CCE COE" "COW COE"; do
  # Split the pair into P1 and P2
  read P1 P2 <<<$(echo $pair)

  # Loop through each iteration
  for IT in {1..5}; do
    # Submit the job with sbatch
    sbatch -J "CROSSCOAL_${P1}_${P2}_${IT}" ~/merondun/cuculus_migration/msmc/4.Crosscoalescent_Iterative.sh ${P1} ${P2} ${IT}
  done
done

```

### Plot MSMC2

```bash
.libPaths('~/mambaforge/envs/R/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/crosscoal/output')
library(tidyverse)
library(viridis)
library(meRo)
library(RColorBrewer)
library(ggpubr)

#for binding later
msmcdat = NULL
files = list.files('.',pattern='.*msmc_ALLHAPS.CROSSCOAL_final.txt')

for (f in files) { 
  cat('Working on file: ',f,'\n')
  t = read_tsv(f)
  string = gsub('_msmc.*','',f)
  split = strsplit(string,'_')[[1]]
  g1 = split[1]
  g2 = split[2]
  run = split[3]
  t$P1 = g1; t$P2 = g2; t$Iteration = run; t$group = paste0(g1,'__',g2)
  msmcdat = rbind(msmcdat,t)
} 

mu = 1.01e-08
gen = 2.74
cols = brewer.pal(12,'Paired')[c(1,2,5,6)]

#label for y axis 
millions_label <- function(x) {
  return(scales::label_number(suffix = "K", scale = 1e-3)(x))
}

#plot each
msmc_each = msmcdat %>% pivot_longer(!c(Iteration,time_index,left_time_boundary,right_time_boundary,P1, P2, group),names_to = 'pop',values_to = 'lambda') %>% 
  mutate(
    #time = left_time_boundary/mu*gen,
    time = right_time_boundary/mu*gen,
    Ne = (1/lambda)/(2*mu)
  )

msmc_each = msmc_each %>% separate(pop,into=c('d1','ID')) %>% 
  mutate(popID = ifelse(ID == '00',P1,
                        ifelse(ID == '11',P2,'Crosspopulation')),
         popID = gsub('D','G',popID)) %>% 
  dplyr::select(-d1,-ID,-P1,-P2)

#remove time boundaries where the estimate for a single population varies high
unstable = msmc_each %>% group_by(popID,time_index) %>% 
  filter(popID != 'Crosspopulation') %>% 
  summarize(mean=mean(Ne),sd=sd(Ne))
unstp = unstable  %>% 
  ggplot(aes(x=time_index,y=sd,col=popID))+
  geom_point(size=0.75)+
  geom_line()+
  scale_color_manual(values=cols)+
  theme_bw()
unstp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/MSMC_UnrelatedBreeding-Instability_2024MAY22.pdf',height=2,width=4)  
unstp
dev.off()

#inspect Ne for outliers
msmc_each %>% filter(popID != 'Crosspopulation' & 
                       time_index >= 3 & 
                       time < 1.5e6) %>% 
  group_by(group, popID) %>%
  ungroup() %>% 
  ggplot(aes(x=popID,fill=group,y=Ne))+geom_boxplot()+theme_bw()

msmc_each %>% count(group)
population_dat = msmc_each %>% filter(popID != 'Crosspopulation' &  
                                        time_index >= 3  & 
                                        time < 1.5e6)
population_dat %>% group_by(group) %>% slice_min(time)
population_dat %>% group_by(group) %>% arrange(desc(Ne))

#harmonic mean
harmonic_mean <- function(x) {
  n <- length(x)
  hm <- n / sum(1 / x)
  return(hm)
}

meanne = population_dat %>% filter(popID != 'Crosspopulation') %>% group_by(popID) %>% summarize(harmNE = harmonic_mean(Ne)) %>% ungroup %>% summarize(lab = paste0(popID,': ',round(harmNE/1000,0),'K',collapse='\n'))
meanne

#plot each 
#ceiling = 2e6
xp = population_dat %>%
  #filter(popID != 'Crosspopulation') %>% 
  #  mutate(Ne = pmin(ceiling,Ne)) %>% 
  ggplot(aes(x = time, y = Ne, col = popID,lty=Iteration)) +
  geom_step(lwd = 0.5, show.legend = TRUE) +
  geom_text(data=meanne,aes(x=Inf,y=Inf,label=lab),vjust=1.5,hjust=1.5,size=1,inherit.aes=FALSE)+
  scale_x_log10(breaks = c(1e4,5e4,1e5,1e6), labels = c('10Ka','50Ka','100Ka','1Ma')) + # Modify breaks as needed
  scale_y_log10(breaks = c(5e4,1e5,2.5e5,5e5,1e6,2.5e6,5e6), labels = c('50K','100K','250K','500K','1M','2.5M','5M')) +  # Modify breaks as needed
  coord_cartesian(xlim=c(1e4,1.5e6),ylim=c(4e4,5e6))+
  theme(strip.text.y = element_text(angle = 0)) +
  scale_color_manual(values=cols)+
  ylab('Ne')+xlab('Time')+
  theme_test(base_size=7)+theme(legend.position='top')
xp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/MSMC_UnrelatedBreeding-NeIterations_2024MAY22.pdf',height=2,width=2)  
xp
dev.off()

#average Ne of each population
population_dat %>% filter(popID != 'Crosspopulation') %>% group_by(popID) %>% sum_stats(Ne)
population_dat %>% filter(popID != 'Crosspopulation') %>% group_by(popID) %>% slice_min(time)
# time_index left_time_boundary right_time_boundary Iteration group    lambda   time       Ne popID
# <dbl>              <dbl>               <dbl> <chr>     <chr>     <dbl>  <dbl>    <dbl> <chr>
#   1          3          0.0000294           0.0000497 5         CCW__CCE  138.  13474.  359867. CCE  
# 2          3          0.0000294           0.0000497 5         CCW__CCE   53.1 13474.  931902. CCW  
# 3          3          0.0000325           0.0000549 1         CCE__COE   33.4 14898. 1484380. COE  
# 4          3          0.0000326           0.0000551 5         CCE__COW   42.8 14938. 1156812. COW 

#or using min/max across runs
popdat_all = population_dat %>%  
  group_by(popID,time_index) %>%
  summarize(time = mean(time),
            minNe = min(Ne),
            maxNe = max(Ne),
            meanNe = mean(Ne))

#plot ribbons 
ribplot = popdat_all %>%
  ggplot(aes(x = time, ymin = minNe, ymax=maxNe, fill = popID)) +
  annotate(geom='rect',xmin=2e4,xmax=2.6e4,ymin=0,ymax=Inf,fill='grey80',alpha=0.5)+
  geom_ribbon(alpha=0.5)+
  #scale_x_log10(breaks = c(1.5e4,2.5e4,5e4,1e5,1e6), labels = c('15Ka','25Ka','50Ka','100Ka','1Ma')) +
  #coord_cartesian(xlim=c(1.5e4,1.5e6),ylim=c(4e4,2.5e6))+
  scale_x_log10(breaks = c(1.5e4,2.5e4,5e4,1e5,2.5e5), labels = c('15Ka','25Ka','50Ka','100Ka','250Ka')) +
  scale_y_log10(breaks = c(5e4,1e5,2.5e5,5e5,1e6,2.5e6,5e6), labels = c('50K','100K','250K','500K','1M','2.5M','5M')) +  # Modify breaks as needed
  coord_cartesian(xlim=c(1e4,1.5e6),ylim=c(4e4,5e6))+
  scale_fill_manual(values=cols)+
  xlab('Time')+ylab('Min Max Ne')+
  theme_test(base_size=7)+
  theme(legend.position='top',
        legend.text=element_text(size=6),
        legend.key.size=unit(0.1,"cm"))

ribplot

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/MSMC_UnrelatedBreeding-NeError250Ka_2024MAY22.pdf',height=2.5,width=2.5)  
ribplot
dev.off()

#plot cross coalesence rate, and add interspecific vs intraspecific comparison 
msmc = msmcdat %>% 
  mutate(time = left_time_boundary/mu*gen,
         crossrate = 2*lambda_01 / (lambda_00 + lambda_11)) %>% 
  mutate(Comparison = ifelse(grepl('CC',P1) & grepl('CC',P2),'C. canorus',
                             ifelse(grepl('CO',P1) & grepl('CO',P2),'C. optatus',
                                    'Interspecific'))) %>% 
  filter(
    time_index >= 3 & 
      time < 1.5e6
  ) 

#or using min/max across runs
cross_all = msmc %>%  
  group_by(Comparison,time_index) %>%
  summarize(time = mean(time),
            mincrossrate = min(crossrate),
            maxcrossrate = pmin(1,max(crossrate)),
            meancrossrate = mean(crossrate))

#when did they split?
splits = cross_all %>% group_by(Comparison) %>% filter(mincrossrate > 0.5) %>% slice_min(time)
splits
# Comparison    time_index   time mincrossrate maxcrossrate meancrossrate
# <chr>              <dbl>  <dbl>        <dbl>        <dbl>         <dbl>
#   1 C. canorus             6 34592.        0.613        0.740         0.686
# 2 C. optatus             5 25189.        0.532        0.632         0.584
# 3 Interspecific          7 60014.        0.504        0.686         0.591

#plot ribbons 
crossrib = cross_all %>%
  ggplot(aes(x = time, ymin = mincrossrate, ymax=maxcrossrate, fill = Comparison)) +
  annotate(geom='rect',xmin=2e4,xmax=2.6e4,ymin=0,ymax=Inf,fill='grey80',alpha=0.5)+
  geom_ribbon(alpha=0.5)+
  #geom_vline(data=splits,aes(xintercept=time,col=Comparison))+
  scale_x_log10(breaks = c(1.5e4,2.5e4,5e4,1e5,2.5e5), labels = c('15Ka','25Ka','50Ka','100Ka','250Ka')) +
  coord_cartesian(xlim=c(1.5e4,2.5e5),ylim=c(0,1))+
  scale_y_continuous(labels=function(x) str_pad(signif(x,2),width=5,side='left'))+
  scale_fill_manual(values=c(cols[c(2,4)],'orange'))+
  scale_color_manual(values=c(cols[c(2,4)],'orange'))+
  xlab('Time')+ylab('Min Max CCR')+
  theme_test(base_size=7)+
  theme(legend.position='top',
        legend.text=element_text(size=6),
        legend.key.size=unit(0.1,"cm"))

crossrib

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/MSMC_UnrelatedBreeding-CrossrateError_2024MAY22.pdf',height=2.5,width=2.5)  
ribplot
dev.off()

#import ice core data 
icecore = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Ice_Temperature_Reconstructions_Kawamura_NOAA-6076.txt')
icecore = icecore %>% dplyr::rename(time = TopAge)

#ice core alone
ic = icecore %>% 
  mutate(filler = '1') %>% 
  ggplot(aes(x = time, y = deltaT,col=filler)) +
  geom_line()+
  annotate(geom='rect',xmin=2e4,xmax=2.6e4,ymin=-Inf,ymax=Inf,fill='grey80',alpha=0.5)+
  scale_y_continuous(labels=function(x) str_pad(signif(x,2),width=5,side='left'))+
  scale_x_log10(breaks = c(1.5e4,2.5e4,5e4,1e5,2.5e5), labels = c('15Ka','25Ka','50Ka','100Ka','250Ka')) + 
  coord_cartesian(xlim=c(1.5e4,2.5e5))+
  scale_color_manual(values=c('black'))+
  xlab('Time')+ylab('Delta T (°C)')+
  theme_test(base_size=7)+
  theme(legend.position='top',
        legend.text=element_text(size=6),
        legend.key.size=unit(0.1,"cm"))
ic

#plot MSMC-IM
imdat = NULL
imfiles = list.files('.',pattern='MSMC_IM.estimates.txt')

for (f in imfiles) { 
  cat('Working on file: ',f,'\n')
  t = read_tsv(f)
  string = gsub('_msmc.*','',f)
  split = strsplit(string,'_')[[1]]
  g1 = split[1]
  g2 = split[2]
  run = split[3]
  t$P1 = g1; t$P2 = g2; t$Iteration = run
  imdat = rbind(imdat,t)
} 

mu = 1.01e-08
gen = 2.74

#add time, drop the first 3 time points as before 
msmcim = imdat %>% 
  group_by(P1, P2, Iteration) %>%
  arrange(left_time_boundary) %>%
  mutate(time_index = row_number(),
         time = left_time_boundary*gen) %>% 
  slice(-c(1:3)) %>%
  ungroup()

#assign group 
msmcim = msmcim %>% 
  mutate(Comparison = ifelse(grepl('CC',P1) & grepl('CC',P2),'C. canorus',
                             ifelse(grepl('CO',P1) & grepl('CO',P2),'C. optatus',
                                    'Interspecific')))

#average by group
mig_dat = msmcim %>% ungroup %>% 
  group_by(Comparison,time_index) %>%
  mutate(time = mean(time)) %>% ungroup %>% 
  group_by(Comparison,time) %>% 
  sum_stats(m)

mp = mig_dat %>%  
  ggplot(aes(x=time,ymin=conf_low,ymax=conf_high,y=mean,col=Comparison,fill=Comparison))+
  geom_line(lwd=0.5)+
  #geom_step(lwd = 0.5, show.legend = TRUE) +
  annotate(geom='rect',xmin=2e4,xmax=2.6e4,ymin=0,ymax=Inf,fill='grey80',alpha=0.5)+
  geom_ribbon(alpha=0.35,col=NA)+
  scale_color_manual(values=c(cols[c(2,4)],'orange'))+
  scale_fill_manual(values=c(cols[c(2,4)],'orange'))+
  scale_x_log10(breaks = c(1.5e4,2.5e4,5e4,1e5,2.5e5), labels = c('15Ka','25Ka','50Ka','100Ka','250Ka')) + 
  coord_cartesian(xlim=c(1.5e4,2.5e5))+
  ylab('Migration Rate')+xlab('Time (Years)')+
  theme_test(base_size=7)+
  theme(legend.position='top',
        legend.text=element_text(size=6),
        legend.key.size=unit(0.1,"cm"))
mp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/MSMC_MigrationRate_2024MAY22.pdf',height=2.5,width=2.5)  
mp
dev.off()


all_plots = ggarrange(ribplot,crossrib,mp,ic,nrow=4)
all_plots

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/MSMC_StackedTimePlots_2024MAY22.pdf',height=6,width=3)  
all_plots
dev.off()

```


