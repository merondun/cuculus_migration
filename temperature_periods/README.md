# Temperature through time

Using the Kawamura ice core data from [here](https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/domefuji/df-tsite-340ka-dfo2006.txt) 

Additional data on the dataset [here](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=noaa-icecore-6076). 

## Assign Two Hot & Cold Periods

Quantitatively assign a hot and cold period, one recent ( < 5e4) and one ancient ( 1.5e5 > x > 1e5 ) 

For example: identifying recent hot times:

```
# Hot times, recent 
hot_young <- icecore %>%  filter(time < 5e4) %>% slice_max(deltaT,n=20) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

  time_low time_high  temp duration
     <dbl>     <dbl> <dbl>    <dbl>
1      500     11500  0.17    11000
```

But, we want our time periods to be a similar duration. I find the midpoitn of the cold periods, and then add +/- 5500 so that all periods are the same duration.

w | time_high | temp | duration | period  | adj_low | adj_high | adj_duration |
| -------- | --------- | ---- | -------- | ------- | ------- | -------- | ------------ |
| 500      | 11500     | hot  | 11000    | recent  | 500     | 11500    | 11000        |
| 18000    | 39500     | cold | 21500    | recent  | 23250   | 34250    | 11000        |
| 121000   | 132000    | hot  | 11000    | ancient | 121000  | 132000   | 11000        |
| 139500   | 143000    | cold | 3500     | ancient | 135750  | 146750   | 11000        |


Which gives us:

![Periods](figures/WarmCold_Periods.png)


Script:

```
.libPaths('~/mambaforge/envs/R/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/crosscoal/output')
library(tidyverse)

#import ice core data 
icecore = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Ice_Temperature_Reconstructions_Kawamura_NOAA-6076.txt')
icecore = icecore %>% dplyr::rename(time = TopAge)

# Hot times, recent 
hot_young <- icecore %>%  filter(time < 5e4) %>% slice_max(deltaT,n=20) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Cold times, recent
cold_young <- icecore %>%  filter(time < 5e4) %>% slice_min(deltaT,n=20) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Hot times, ancient
hot_old <- icecore %>%  filter(time > 1e5 & time < 1.5e5) %>% slice_max(deltaT,n=20) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

# Cold times, ancient
cold_old <- icecore %>%  filter(time > 1e5 & time < 1.5e5) %>% slice_min(deltaT,n=3) %>% 
  summarize(time_low = min(time), time_high = max(time), temp = mean(deltaT), duration = time_high - time_low)

periods <- rbind(hot_young, cold_young, hot_old, cold_old) %>% 
  mutate(temp = c('hot','cold','hot','cold'),
         period = c('recent','recent','ancient','ancient'))
periods

# Hot times are 11 Ka, but our cold times are 21.5 Ka and 19 Ka, so just take midpoint of those cold dates
periods <- periods %>% 
  mutate(
    #Calculate midpoint (duration / 2 + low point, and then add the interval on either side 
    adj_low = ifelse(temp == 'cold', ((duration / 2 ) + time_low) - 5500, time_low),
    adj_high = ifelse(temp == 'cold', ((duration / 2 ) + time_low) + 5500, time_high),
    adj_duration = adj_high - adj_low
  )

# Kawamura NOAA 6076 Ice core inferred temperature
ice_plot <- icecore %>% 
  ggplot(aes(x = time, y = deltaT)) +
  geom_line()+
  #annotate(geom='rect',xmin=2e4,xmax=2.6e4,ymin=-Inf,ymax=Inf,fill='grey80',alpha=0.5)+
  geom_text(data=periods,aes(x=adj_low+5500,y=4,label=temp),col='black')+
  coord_cartesian(xlim=c(0,1.75e5))+
  geom_rect(data=periods,aes(xmin=adj_low,xmax=adj_high,ymin=-Inf,ymax=Inf,fill=temp),
            alpha=0.5,inherit.aes=FALSE)+
  scale_fill_manual(values=c('cyan3','salmon2'))+
  xlab('Time')+ylab('Delta T (°C)')+
  theme_test(base_size=7)
ice_plot

png('~/merondun/cuculus_migration/figures/WarmCold_Periods.png',quote=F,sep='\t',row.names=F)
ice_plot
dev.off()
```
