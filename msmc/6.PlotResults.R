.libPaths('~/mambaforge/envs/R/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/crosscoal/output')
library(tidyverse)
library(viridis)
library(meRo)
library(RColorBrewer)
library(gmodels)

#for binding later
msmcdat = NULL
files = list.files('.',pattern='ALLHAPS.CROSSCOAL_final.txt')

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
  scale_color_viridis(discrete=TRUE)+
  theme_bw(base_size=8)
unstp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_3pop_All-Instability.pdf',height=2,width=3)  
unstp
dev.off()

#inspect Ne for outliers
msmc_each %>% 
  filter(popID != 'Crosspopulation' & time_index >= 3) %>% 
  group_by(group, popID) %>%
  ungroup() %>% 
  ggplot(aes(x=popID,fill=group,y=Ne))+geom_boxplot()+theme_bw()

msmc_each %>% count(group)

# Filter for final plotting, removing the highly unstable time periods (<3)
population_dat = msmc_each %>% filter(popID != 'Crosspopulation' & time_index >= 3)

cols <- data.frame(popID = c('CCW','CCE','CO','COW'),
                   col = c('#33a02c','#b2df8a','gold','#fe9024'))

#plot each 
xp = population_dat %>%
  #filter(popID != 'COW') %>%  
  ggplot(aes(x = time, y = Ne, col = popID,lty=Iteration)) +
  geom_step(lwd = 0.25, show.legend = TRUE) +
  scale_x_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,1e6), labels = c('10KY','100KY','1MY')) + # Modify breaks as needed
  scale_y_log10(limits=c(1e4,3e6),breaks = c(1e4,5e4,1e5,2.5e5,5e5,1e6), labels = millions_label) +  # Modify breaks as needed
  theme(strip.text.y = element_text(angle = 0)) +
  scale_color_manual(values=cols$col,breaks=cols$popID)+
  ylab('Ne')+xlab('Time')+
  theme_test(base_size=7)+theme(legend.position='top')
xp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_3pop_All_NeLines.pdf',height=2,width=2)  
xp
dev.off()

#average Ne of each population
population_dat %>% filter(popID != 'Crosspopulation') %>% group_by(popID) %>% sum_stats(Ne)
population_dat %>% filter(popID != 'Crosspopulation') %>% group_by(popID) %>% slice_min(time)

# popID    mean      sd     se  median    iqr conf_low conf_high
# <chr>   <dbl>   <dbl>  <dbl>   <dbl>  <dbl>    <dbl>     <dbl>
#   1 CCE   105574.  77034.  3596.  90950. 69297.   98508.   112640.
# 2 CCW   190236. 250428. 11689. 124133. 74881.  167265.   213206.
# 3 CO    243492. 363291. 20768. 156906. 93296.  202625.   284359.
# 4 COW   219910. 281023. 16065. 157642. 91971.  188298.   251523.
# # 
# time_index left_time_boundary right_time_boundary Iteration group    lambda   time       Ne popID
# <dbl>              <dbl>               <dbl> <chr>     <chr>     <dbl>  <dbl>    <dbl> <chr>
#   1          3          0.0000292           0.0000494 1         CCW__CCE  187.  13395.  264610. CCE  
# 2          3          0.0000292           0.0000494 1         CCW__CCE   44.0 13395. 1126166. CCW  
# 3          3          0.0000322           0.0000544 1         CCE__CO    29.1 14761. 1698615. CO   
# 4          3          0.0000325           0.0000549 4         CCE__COW   35.1 14898. 1411688. COW  

#or using min/max across runs
popdat_all = msmc_each %>% filter(popID != 'Crosspopulation' & time_index >= 3 & time != Inf) %>% 
  group_by(popID,time_index) %>%
  summarize(time = min(time),
            minNe = min(Ne),
            maxNe = max(Ne),
            meanNe = mean(Ne))

#plot ribbons 
ribplot = popdat_all %>%
  ggplot(aes(x = time, ymin = minNe, ymax=maxNe, fill = popID)) +
  geom_ribbon(alpha=0.5)+
  scale_x_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,1e6), labels = c('10KY','100KY','1MY')) + # Modify breaks as needed
  scale_y_log10(limits=c(1e4,3e6),breaks = c(1e4,5e4,1e5,2.5e5,5e5,1e6), labels = millions_label) +  # Modify breaks as needed
  theme(strip.text.y = element_text(angle = 0)) +
  scale_fill_manual(values=cols$col,breaks=cols$popID)+
  ylab('Min and Max Ne')+xlab('Time')+
  theme_test(base_size=7)
ribplot

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_3pop-All_Ribbons.pdf',height=1.5,width=2.5)  
ribplot
dev.off()


# #interplotated to smooth ribbons 
# n_interpolated_points = 10000
# popdat_all_interpolated <- popdat_all %>%
#   group_by(popID) %>% 
#   mutate(
#     time_new = list(seq(min(time), max(time), length.out = n_interpolated_points)),
#     minNe_new = list(spline(time, minNe, n = n_interpolated_points)$y),
#     maxNe_new = list(spline(time, maxNe, n = n_interpolated_points)$y)
#   ) %>%
#   ungroup() %>%
#   select(-time, -minNe, -maxNe) %>%
#   unnest(c(time_new, minNe_new, maxNe_new)) %>%
#   rename(time = time_new, minNe = minNe_new, maxNe = maxNe_new)
# 
# ribplot = popdat_all_interpolated %>%
#   ggplot(aes(x = time, ymin = minNe, ymax=maxNe, fill = popID)) +
#   #geom_line(data=popdat_all,aes(x=time,y=meanNe))+
#   geom_ribbon(alpha=0.8)+
#   scale_x_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,1e6), labels = c('10KY','100KY','1MY')) + # Modify breaks as needed
#   scale_y_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,2.5e5,5e5,1e6), labels = millions_label) +  # Modify breaks as needed
#   theme(strip.text.y = element_text(angle = 0)) +
#   scale_fill_viridis(discrete=TRUE)+
#   ylab('Min and Max Ne')+xlab('Time')+
#   theme_test(base_size=7)+
#   theme(legend.position='top')
# ribplot
# 
# pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/figures/MSMC_FourPopulations_Introgression-Error_2024APR1.pdf',height=2.5,width=2.5)  
# ribplot
# dev.off()

#plot cross coalesence rate, and add interspecific vs intraspecific comparison 
msmc = msmcdat %>% 
  #filter(time_index >= 3) %>% 
  #filter(P1 != 'CO' & P2 != 'CO') %>% 
  mutate(time = left_time_boundary/mu*gen,
         crossrate = pmin(1,2*lambda_01 / (lambda_00 + lambda_11))) %>% 
  mutate(Comparison = ifelse(grepl('CC',P1) & grepl('CC',P2),'C. canorus',
                             ifelse(grepl('CO',P1) & grepl('CO',P2),'C. optatus',
                                    'Interspecific')))

#import ice core data 
icecore = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Ice_Temperature_Reconstructions_Kawamura_NOAA-6076.txt')
icecore = icecore %>% dplyr::rename(time = TopAge)

#scale ice core according to the cross rate
crossrate_range <- range(msmc$crossrate, na.rm = TRUE)
deltaT_range <- range(icecore$deltaT, na.rm = TRUE)
scale_factor <- diff(crossrate_range) / diff(deltaT_range)
icecore$scaled_deltaT <- (icecore$deltaT - mean(deltaT_range)) * scale_factor + mean(crossrate_range)
msmc <- msmc %>% mutate(Opt = ifelse(grepl('COW$',P1) | grepl('COW$',P2),'C. optatus western (n=10)','Merged C. optatus (n=20)'))

climates <- data.frame(
  xmin = c(0, 24e3, 57e3, 122e3),        
  xmax = c(10e3, 34e3, 67e3, 132e3),   
  ymin = -Inf,                           
  ymax = Inf,                            
  fill = c("#c3e6ee","#ffbfbf","#ffbfbf","#c3e6ee"))

#  fill = c("#ffbfbf","#c3e6ee","#c3e6ee", "#ffbfbf"))

ap = msmc %>% 
  #filter(Iteration == 1) %>% 
  ggplot(aes(x = time, y = crossrate,Group=interaction(P1,P2),col=interaction(P1,P2),lty=Iteration)) +
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.1, inherit.aes=FALSE) +
  geom_step(lwd = 0.25, show.legend = TRUE) +
  scale_x_log10(breaks = c(1e4,2.5e4,5e4,1e5,1.5e5,2.5e5), labels = c('10 Ka','25 Ka','50 Ka','100 Ka','150 Ka','250 Ka')) + 
  coord_cartesian(ylim=c(0,1),xlim=c(2e3,3.5e5))+
  scale_y_continuous(
    name = "Crossrate",
    sec.axis = sec_axis(~ . / scale_factor + mean(deltaT_range), name = "Delta T (°C)")
  ) +
  theme(strip.text.y = element_text(angle = 0)) +
  geom_line(data = icecore, aes(x = time, y = scaled_deltaT, group = 1), lwd=0.25,lty=2,inherit.aes=FALSE,color = "black")+
  scale_color_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c'))+
  theme_test(base_size=7)+
  facet_grid(Opt~.)+
  theme(legend.position='top')
ap

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_CrosscoalLines-3pop-All_Migration_IceCore.pdf',height=3,width=3)  
ap
dev.off()

## Rib plot for crosscoal
#or using min/max across runs
popdat_ccr = msmc %>% 
  group_by(time_index,P1,P2,Comparison,Opt) %>%
  mutate(time = min(time)) %>% ungroup %>% 
  group_by(time,P1,P2,Comparison,Opt) %>% 
  sum_stats(crossrate)

#plot ribbons 
ribplot_ccr = popdat_ccr %>%
  ggplot(aes(x = time, ymin = conf_low, Group=interaction(P1,P2), ymax=conf_high, fill = interaction(P1,P2))) +
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.3, inherit.aes=FALSE) +
  geom_ribbon(alpha=0.5)+
  scale_x_log10(breaks = c(1e4,2.5e4,5e4,1e5,1.5e5,2.5e5), labels = c('10 Ka','25 Ka','50 Ka','100 Ka','150 Ka','250 Ka')) + 
  coord_cartesian(ylim=c(0,1),xlim=c(2e3,3.5e5))+
  scale_y_continuous(
    name = "Crossrate",
    sec.axis = sec_axis(~ . / scale_factor + mean(deltaT_range), name = "Delta T (°C)")
  ) +
  scale_fill_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c',"#c3e6ee","#ffbfbf"))+
  theme(strip.text.y = element_text(angle = 0)) +
  geom_line(data = icecore, aes(x = time, y = scaled_deltaT, group = 1),lwd=0.25, lty=2,inherit.aes=FALSE,color = "black")+
  scale_color_manual(values=brewer.pal(3,'Set2'))+
  theme_test(base_size=7)+
  facet_grid(Opt~.)+
  theme(legend.position='top')
ribplot_ccr

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_CrosscoalRibbons-3pop-All_Migration_IceCore.pdf',height=2.5,width=3)  
ribplot_ccr
dev.off()


#plot MSMC-IM
imdat = NULL
imfiles = list.files('.',pattern='MSMC_IM.estimates.txt')
#imfiles <- imfiles[!grepl('_CO_',imfiles)]

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
imdat <- imdat %>% mutate(Opt = ifelse(grepl('COW$',P1) | grepl('COW$',P2),'C. optatus western (n=10)','Merged C. optatus (n=20)'))
msmcim = imdat %>% 
  group_by(P1, P2, Iteration,Opt) %>%
  arrange(left_time_boundary) %>%
  mutate(time_index = row_number(),
         time = left_time_boundary*gen) %>%  #note that mu is already integrated into MSMC-IM! 
  #slice(-c(1:4)) %>%
  ungroup() %>% 
  mutate(Comparison = ifelse(grepl('CC',P1) & grepl('CC',P2),'C. canorus',
                             ifelse(grepl('CO',P1) & grepl('CO',P2),'C. optatus',
                                    'Interspecific')))

mp = msmcim %>%  
  ggplot(aes(x=time,y=m,group=interaction(P1,P2,Iteration),lty=Iteration,col=interaction(P1,P2)))+
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.01, inherit.aes=FALSE) +
  geom_step(show.legend = TRUE) +
  scale_color_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c'))+
  #facet_wrap(~Opt,scales='free')+
  scale_x_log10(limits=c(5e3,1e5),breaks = c(5e3,1e4,2e4,5e4,1e5), labels = c('5 Ka','10 Ka','20 Ka','50Ka','100 Ka'))+
  ylab('Migration Rate (m)')+xlab('Time (Years)')+
  theme_test(base_size=7)+
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.2, inherit.aes=FALSE) +
  facet_grid(Opt~.)+
  theme(legend.position='top')
mp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMCIM_3pop-All_AllLines.pdf',height=2.5,width=3)  
mp
dev.off()


#average by group
mig_dat = msmcim %>% ungroup %>% 
  group_by(P1,P2,Comparison,time_index,Opt) %>%
  mutate(time = min(time)) %>% ungroup %>% 
  group_by(time,P1,P2,Comparison,Opt) %>% 
  sum_stats(m) %>% 
  ungroup %>% 
  arrange(time)

#and just plot 1 for final figure 
mp2 = mig_dat %>%  
  ggplot(aes(x=time,ymin=conf_low,ymax=conf_high,group=interaction(P1,P2),fill=interaction(P1,P2)))+
  geom_ribbon(alpha=0.5)+
  scale_x_log10(limits=c(3000,1e5),breaks = c(5e3,1e4,2e4,5e4,1e5), labels = c('5 Ka','10 Ka','20 Ka','50Ka','100 Ka'))+
  ylab('Migration Rate (m)')+xlab('Time (Years)')+
  theme_test(base_size=7)+
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.3, inherit.aes=FALSE) +
  facet_grid(Opt~.)+
  scale_fill_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c',"#ffbfbf","#c3e6ee"))+
  theme(legend.position='top')
mp2

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMCIM_3pops_All-Ribbons.pdf',height=2.5,width=2.5)  
mp2
dev.off()

