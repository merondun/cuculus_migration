.libPaths('~/mambaforge/envs/R/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/crosscoal/output')
library(tidyverse)

#import ice core data 
icecore = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Ice_Temperature_Reconstructions_Kawamura_NOAA-6076.txt')
icecore = icecore %>% dplyr::rename(time = TopAge)

# Hot times, recent 
hot_0 <- icecore %>%  filter(time < 5e4) %>% slice_max(deltaT,n=20) %>% 
  summarize(time_low = 0, time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Cold times, recent
cold_0 <- icecore %>%  filter(time < 5e4) %>% slice_min(deltaT,n=20) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Hot times, mid
hot_1 <- icecore %>%  filter(time < 1.75e5 & time > 5e4) %>% slice_max(deltaT,n=20) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Cold times, mid
cold_1 <- icecore %>%  filter(time < 1e5 & time > 5e4) %>% slice_min(deltaT,n=20) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

periods <- rbind(hot_0, cold_0, hot_1, cold_1) %>% 
  mutate(stage = c('hot','cold','hot','cold'),
         period = c('t0','t0','t1','t1'))
periods

# Adjusting durations to be exactly 10000
adjusted_periods <- periods %>%
  mutate(
    mid_point = (time_low + time_high) / 2,  #Calculate midpoint
    time_low = if_else(duration > 10000, mid_point - 5000, time_low), 
    time_high = if_else(duration > 10000, mid_point + 5000, time_high),  
    duration = 10000  # Set duration to 10000
  ) %>%
  select(-mid_point)  # Remove the mid_point column

# Output the adjusted dataframe
adjusted_periods

#Modify recent time so that it occurs from contemporary time 
adjusted_contemp <- adjusted_periods %>% mutate(
  time_low = ifelse(stage == 'hot' & period == 't0',0,time_low),
  time_high = ifelse(stage == 'hot' & period == 't0',10000,time_high)) %>% 
  select(-temp) %>% #recalculate temps within this 10KY
  rowwise() %>%
  mutate(
    temp = mean(icecore$deltaT[icecore$time >= time_low & icecore$time <= time_high], na.rm = TRUE)
  )

adjusted_contemp

# Kawamura NOAA 6076 Ice core inferred temperature
ice_plot <- icecore %>% 
  ggplot(aes(x = time, y = deltaT)) +
  geom_line()+
  coord_cartesian(xlim=c(0,1.5e5))+
  geom_rect(data=adjusted_contemp,aes(xmin=time_low,xmax=time_high,ymin=-Inf,ymax=Inf,fill=stage),
            alpha=0.5,inherit.aes=FALSE)+
  scale_fill_manual(values=c('cyan3','salmon2'))+
  geom_label(data=adjusted_contemp,size=1.5,aes(x=time_low+5500,y=4,
                                                label=paste0(stage,': ',period,'  ',round(temp,2),'°C\n',time_low/1000,' - ',time_high/1000,'Ka')),col='black')+
  xlab('Time')+ylab('Delta T (°C)')+
  theme_test(base_size=10)
ice_plot

png('~/merondun/cuculus_migration/figures/WarmCold_Periods.png',units='in',res=300,height=2,width=3.5)
pdf('~/merondun/cuculus_migration/figures/WarmCold_Periods.pdf',height=2,width=3.5)
ice_plot
dev.off()
