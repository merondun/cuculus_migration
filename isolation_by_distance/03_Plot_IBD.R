#### Plot Distance ~ FST
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/isolation_by_distance')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(broom)
library(ggpubr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(meRo)
library(RColorBrewer)

#Read in metadata
md = read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')

d = read_tsv('20250122_Pairwise_GeographicDistance_Km.txt')
f = read_tsv('20250122_Pairwise_FST_DistanceGroups_Autosomes.txt')

#prep and merge frames 
names(f) = c('chr','start','end','snps','FST','Group')
f2 = f %>% mutate(FST = pmax(0, pmin(1, FST))) %>% select(!c(start,end,snps))

#merge with geographic distance 
df = left_join(f2,d) 

dfs = df %>% group_by(Group,P1,P2,Distance_km) %>% 
  sum_stats(FST)
dfs %>% ggplot(aes(x=log(Distance_km),y=log(mean)))+
  geom_point()+
  geom_smooth()+
  theme_bw()

dfs <- dfs %>% mutate(Species = ifelse(grepl('CO',Group),'C. optatus','C. canorus'))

# Spearman's rho and cor test within each AvZ level and gather coefficients
cors = dfs %>% group_by(Species) %>%
  summarize(
    rho = cor(mean, Distance_km,method='spearman'),
    p_value = cor.test(mean, Distance_km,method='spearman')$p.value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'),
         signif = case_when(
           padj < 0.05 ~ "*",
           TRUE ~ " (n.s.)"),
         label = paste0('rs: ',round(rho, 2), signif))
cors
# Species      rho p_value  padj signif    label          
# <chr>      <dbl>   <dbl> <dbl> <chr>     <chr>          
#   1 C. canorus 0.772   0     0     "*"       rs: 0.77*      
#   2 C. optatus 0.221   0.427 0.853 " (n.s.)" rs: 0.22 (n.s.)

cols = c('darkorchid2','orange')
dfs$Species = factor(dfs$Species,levels=c('C. canorus','C. optatus'))

#lmm, account for P1 and P2 
model_summaries <- dfs %>%
  group_by(Species) %>%
  do({
    model <- lmer(log(pmax(mean, 0.005)) ~ log(pmax(Distance_km, 0.005)) + (1|P1) + (1|P2), data = .)
    summary_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  }) %>%
  dplyr::bind_rows()

#output fixed effect of distance, bonferroni correction 
model_summaries %>% filter(grepl('Distance',term)) %>% select(-effect,-term) %>% ungroup %>% 
  mutate(padj = p.adjust(p.value,method='bonferroni'),
         signif = ifelse(padj < 0.05,'*','n.s.'))

# Species    estimate std.error statistic    df  p.value conf.low conf.high     padj signif
# <fct>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr> 
#   1 C. canorus    0.743    0.0465     16.0  98.8  3.82e-29    0.651     0.836 7.64e-29 *     
#   2 C. optatus    0.446    0.0691      6.45  5.43 9.77e- 4    0.272     0.619 1.95e- 3 *  


# Create the plot
pp1 = ggplot(dfs, aes(x = log(pmax(Distance_km,0.005)), y = log(pmax(mean,0.005)), color = Species, shape = Species)) +
  geom_point(size=1) + #0.25 for main plot 
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +  #0.5 for main plot 
  geom_text(data = cors[1,], aes(x = -Inf, y = Inf, label = paste0('C. canorus rho: ',round(rho,3))), size = 3,vjust = 2, hjust = -.25) +
  geom_text(data = cors[2,], aes(x = -Inf, y = Inf, label = paste0('C. optatus rho: ',round(rho,3))), size = 3,vjust = 4, hjust = -.25) +
  labs(title = "Relationship between FST and Distance_km",
       x = "log(Geographic Distance (km))",
       y = "log(FST)") +
  scale_color_manual(values=cols)+
  theme_bw(base_size=6) + theme(legend.position = 'none')
pp1


ggsave('~/merondun/cuculus_migration/figures/20250115_FSTvGeoDistanceSpecies.pdf',pp1,height=3,width=3.5,dpi=300)


#### Fst / (1 - FST) and log(dist)

cors_rousset = dfs %>% group_by(Species) %>%
  summarize(
    rho = cor(mean, Distance_km,method='spearman'),
    p_value = cor.test( (mean / (1 - mean)), Distance_km,method='spearman')$p.value) %>% 
  mutate(label = paste0('rs: ',round(rho, 2)))
cors_rousset

#lmm, account for P1 and P2 
model_summaries_rousset <- dfs %>%
  group_by(Species) %>%
  do({
    model <- lmer( (mean / (1 - mean)) ~ log(pmax(Distance_km, 0.005)) + (1|P1) + (1|P2), data = .)
    summary_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  }) %>%
  dplyr::bind_rows()

#output fixed effect of distance, bonferroni correction 
model_summaries_rousset %>% filter(grepl('Distance',term)) %>% select(-effect,-term) %>% ungroup %>% 
  mutate(padj = p.adjust(p.value,method='bonferroni'),
         signif = ifelse(padj < 0.05,'*','n.s.'))

# Species    estimate std.error statistic    df  p.value conf.low conf.high     padj signif
# <fct>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr> 
#   1 C. canorus  0.0310   0.00220      14.1  89.9  1.85e-24  0.0266    0.0354  3.70e-24 *     
#   2 C. optatus  0.00737  0.000740      9.96  5.14 1.49e- 4  0.00548   0.00926 2.98e- 4 *   


# Create the plot
pp2 = ggplot(dfs, aes(x = log(pmax(Distance_km,0.005)), y = (mean / (1 - mean)), color = Species, shape = Species)) +
  geom_point(size=1) + #0.25 for main plot 
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +  #0.5 for main plot 
  geom_text(data = cors[1,], aes(x = -Inf, y = Inf, label = paste0('C. canorus rho: ',round(rho,3))), size = 3,vjust = 2, hjust = -.25) +
  geom_text(data = cors[2,], aes(x = -Inf, y = Inf, label = paste0('C. optatus rho: ',round(rho,3))), size = 3,vjust = 4, hjust = -.25) +
  labs(title = "Relationship between FST and Distance_km",
       x = "log(Geographic Distance (km))",
       y = "FST / (1 - FST)") +
  scale_color_manual(values=cols)+
  theme_bw(base_size=6) + theme(legend.position = 'none')
pp2


ggsave('~/merondun/cuculus_migration/figures/20250115_FSTvGeoDistanceRoussetSpecies.pdf',pp1,height=3,width=3.5,dpi=300)


