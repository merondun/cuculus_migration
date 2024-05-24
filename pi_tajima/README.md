# Estimating π & Tajima's D

Estimate Tajima's D and π in 25Kb windows, including invariant sites. Calculate with [genomics_general](https://github.com/simonhmartin/genomics_general/tree/master). Note that FST estimated from Simon Martin's scripts is the *Kst* statistic from [this paper](https://doi.org/10.1093/oxfordjournals.molbev.a040703) (e.g. F<sub>ST</sub> = 1 - π<sub>S</sub>  / π<sub>T</sub> ), or pairwise diversity within pops compared to total diversity. 

![Pi Plot](/figures/Pi_TajimasD.png)

Highest diversity in the C. optatus populations, with the lowest in C. canorus east. Lowest Tajima's D in C. canorus west supports strongest expansion. 


| Group | Vawriable   | mean    | sd       | se       | median | iqr    | 95ci_low | 95ci_high |
| ----- | ----------- | ------- | -------- | -------- | ------ | ------ | -------- | --------- |
| CCE   | Tajimas's D | -0.486  | 0.548    | 0.0027   | -0.467 | 0.655  | -0.491   | -0.48     |
| CCW   | Tajimas's D | -0.762  | 0.419    | 0.00206  | -0.763 | 0.493  | -0.766   | -0.758    |
| COE   | Tajimas's D | -0.556  | 0.432    | 0.00213  | -0.568 | 0.51   | -0.56    | -0.552    |
| COW   | Tajimas's D | -0.545  | 0.441    | 0.00217  | -0.563 | 0.529  | -0.55    | -0.541    |
| CCE   | π           | 0.00123 | 0.000528 | 2.6E-06  | 0.0012 | 0.0008 | 0.00122  | 0.00123   |
| CCW   | π           | 0.00146 | 0.000564 | 2.77E-06 | 0.0015 | 0.0008 | 0.00146  | 0.00147   |
| COE   | π           | 0.00163 | 0.000592 | 2.91E-06 | 0.0017 | 0.0008 | 0.00162  | 0.00163   |
| COW   | π           | 0.00164 | 0.000604 | 2.97E-06 | 0.0017 | 0.0008 | 0.00163  | 0.00164   |



```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00
#SBATCH --output=/dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/pi_tajima/Slurms/slurm-%j.out

#mamba activate snps
# Submit as: for CHR in $(cat /dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/Chromosomes.list); do sbatch -J PI_${CHR} 1_EstimatePiTajima.sh ${CHR}; done  

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
```

Simply plot in R:

```R
#### pi & tajima's D estimates
setwd('~/merondun/cuculus_migration/pi_tajima/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(meRo)

#read in data 
cols <- brewer.pal(12,'Paired')[c(2,1,6,5)]
df <- read_tsv('Pi_Tajima_2024MAY23.txt.gz',col_names = T)

# Select columns and gather them into a 'long' format
pi <- df %>%
  select(scaffold, start, end, starts_with("pi_")) %>%
  gather(key = "Group", value = "π", -scaffold, -start, -end) %>% 
  group_by(Group) %>% 
  mutate(Group = gsub('pi_','',Group))

taj <- df %>%
  select(scaffold, start, end, starts_with("TajD")) %>%
  gather(key = "Group", value = "Tajimas's D", -scaffold, -start, -end,) %>% 
  group_by(Group) %>% 
  mutate(Group = gsub('TajD_','',Group))

# Mergee data
tapi <- full_join(pi,taj)

# Only show scientific notation for pi because it's tiny
format_numbers <- function(x, threshold = 0.01) {
  sapply(x, function(num) {
    if (abs(num) < threshold) {
      format(num, scientific = TRUE)
    } else {
      format(num, scientific = FALSE)
    }
  })
}

# Plot 95% CIs
sums <- tapi %>%
  pivot_longer(c(π,"Tajimas's D")) %>% 
  group_by(Group,name) %>%
  sum_stats(value) %>%
  mutate(label = paste0(format_numbers(signif(median,3)),' ± ',format_numbers(signif(iqr,3))))
tapi$Group <- factor(tapi$Group,levels=c('CCW','CCE','COW','COE'))

# Plot
tapi_plot <- tapi %>% 
  pivot_longer(c(π,"Tajimas's D")) %>% 
  ggplot(aes(x=Group,fill=Group,y=value))+
  ggdist::stat_halfeye(width = .3,.width = 0, trim=FALSE, justification = 0, 
                       point_colour = NA,normalize='groups')+
  geom_boxplot(width=0.25,position=position_nudge(x=-0.25),outlier.alpha = 0.2)+
  geom_text(data=sums,aes(x=Group,y=Inf,label=label),inherit.aes=FALSE,vjust=1.5,size=1.5)+
  facet_wrap(name~.,scales='free')+
  ylab('Value')+
  scale_fill_manual(values=cols)+
  theme_bw(base_size=8)
tapi_plot

png('../figures/Pi_TajimasD.png',units='in',res=300,height=3,width=6.5)
tapi_plot
dev.off()
```
