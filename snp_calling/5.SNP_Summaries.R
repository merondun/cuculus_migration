#### SNP Stats 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs/stats')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)

d = read_tsv('Filtering_SNP_Counts_2023DEC20.txt',col_names = F)
names(d) = c('chr','Filter','SNPs')
d = d %>% mutate(chr = gsub('chr_','',chr))
d$Filter = factor(d$Filter,levels=c('LD_R2_W50','Polarized','5X_MM1','MQ20','AllSites','Raw'))
chord = d %>% select(chr) %>% unique %>% arrange(chr) %>% mutate_at('chr',as.numeric) %>% arrange(chr) 
d$chr = factor(d$chr,levels=chord$chr)
ds = d %>% group_by(Filter) %>% 
  summarize(SNPs = sum(SNPs))
ds %>% ggplot(aes(x=Filter,y=SNPs,fill=Filter))+
  geom_text(aes(label=paste0(signif(SNPs/1e6,2),'M')),vjust=-1)+
  geom_bar(stat='identity')+
  theme_bw()
df = d %>% filter(!grepl('Raw|AllSites',Filter))
df$Filter = factor(df$Filter,levels=c('MQ20','5X_MM1','Polarized','LD_R2_W50'))

#count proportions
df = df %>%
  group_by(chr) %>%
  mutate(total_SNPs = sum(SNPs[Filter == "MQ20"]),
         proportion = SNPs / total_SNPs,
         SNPs_million = signif(SNPs / 1e4, 2)) %>%
  ungroup()

#bars
ggplot(df, aes(x = chr, y = proportion, fill = Filter)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(aes(label = paste0(SNPs_million,'K'),y=proportion), 
            position = position_stack(vjust=-0.05), 
            color = "black",size=2) +
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Chromosome", y = "Proportion of SNPs", fill = "Filter") +
  theme_minimal()

#total just by filter 
df %>% group_by(Filter) %>% summarize(SNPs_M = round(sum(SNPs_million)/1000,2)) %>% 
  ggplot(aes(x = Filter, y = SNPs_M,fill=Filter)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(aes(label = paste0(SNPs_M,'M'),y=Inf,vjust=1), 
            color = "black",size=3) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Filter", y = "SNPS Retained", fill = "Filter") +
  theme_minimal()
