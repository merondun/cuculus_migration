.libPaths('~/mambaforge/envs/R/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/crosscoal/output')
setwd('G:/My Drive/Research/Migration/temperature')
library(tidyverse)

#import ice core data 
icecore = read_tsv('Ice_Temperature_Reconstructions_Kawamura_NOAA-6076.txt')
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
                                                label=paste0(stage,': ',period,'  ',round(temp,2),'°C/n',time_low/1000,' - ',time_high/1000,'Ka')),col='black')+
  xlab('Time')+ylab('Delta T (°C)')+
  theme_test(base_size=10)
ice_plot

png('~/merondun/cuculus_migration/figures/WarmCold_Periods.png',units='in',res=300,height=2,width=3.5)
pdf('~/merondun/cuculus_migration/figures/WarmCold_Periods.pdf',height=2,width=3.5)
ice_plot
dev.off()

# Correlation with pastclim
install.packages("pastclim")
library(pastclim)
get_available_datasets()
get_vars_for_dataset(dataset = "Beyer2020")
set_data_path(path_to_nc = "G:/My Drive/Research/Migration/temperature/")
download_dataset(dataset = "Beyer2020", bio_variables = c("bio01"),annual=TRUE,monthly=FALSE)
# Get all available time steps for Beyer2020
time_steps <- get_time_bp_steps(dataset = "Beyer2020")

# Extract the global raster series for bio01 across all time steps
global_bio01 <- region_series(
  dataset = "Beyer2020",
  bio_variables = "bio01",
  time_bp = time_steps # Provide explicit vector of time steps
)
# Ensure the SpatRaster contains all layers
global_means <- sapply(1:terra::nlyr(global_bio01$bio01), function(layer_index) {
  layer <- terra::subset(global_bio01$bio01, layer_index) # Get layer by index
  terra::global(layer, fun = "mean", na.rm = TRUE)$mean # Compute mean temperature
})

# Combine results with time steps
results <- data.frame(
  Time_BP = time_steps, # Use the time_steps vector
  Global_Mean_Temperature = global_means
)

combined <- left_join(results %>% mutate(Time_BP = abs(Time_BP)),icecore %>% dplyr::rename(Time_BP = time)) %>% na.omit
temp_cor <- cor(combined$Global_Mean_Temperature,combined$deltaT,method='spearman')
temp_plot <- combined %>% 
  ggplot(aes(x=Global_Mean_Temperature,y=deltaT))+
  geom_point()+xlab('Temperature, Beyer 2020')+ylab('Temperature Change, Kawamura 2007')+
  ggtitle(paste0('Correlation Between Hindcasting & \nDemographic Temperature Datasets (rho = ',round(temp_cor,3),')'))+
  theme_bw()
ggsave(temp_plot,file='Temperature_Correlation_Beyer2020-Kawamura2007.png',dpi=300,height=5,width=5)

# Plot Periods from Beyer2020
combined <- right_join(results %>% mutate(Time_BP = abs(Time_BP)),icecore %>% dplyr::rename(Time_BP = time))
beyer_plot <- combined %>% 
  mutate(deltaT = deltaT + 12) %>% 
  dplyr::rename(Beyer2020 = Global_Mean_Temperature, Kawamura2007 = deltaT) %>% 
  pivot_longer(!Time_BP, names_to = 'Dataset', values_to = 'Temperature') %>% 
  na.omit %>% 
  ggplot(aes(x = Time_BP, y = Temperature, lty = Dataset)) +
  geom_line()+
  coord_cartesian(xlim=c(0,1.5e5))+
  geom_rect(data=adjusted_contemp,aes(xmin=time_low,xmax=time_high,ymin=-Inf,ymax=Inf,fill=stage),
            alpha=0.5,inherit.aes=FALSE)+
  scale_fill_manual(values=c('cyan3','salmon2'))+
  xlab('Time')+ylab('Delta T (°C)')+
  theme_test(base_size=10)
beyer_plot

ggsave(beyer_plot,file='Temperature_Periods_Beyer2020.png',dpi=300,height=4,width=8)

