### Subset N = 80 based on geography. 
setwd('~/merondun/cuculus_migration/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(terra)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(meRo) #devtools::install_github('merondun/meRo')
library(rworldmap)
library(sf) #for spatial plotting of distance vectors to confirm
library(ggspatial) #to add scale bars onto maps
library(ggpubr)
library(lubridate)
library(geosphere)

#read in metadata
md = read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')

#popgen dat 
all = md %>% filter(Analysis_ADMIXTURE_Breeding == 1 ) %>% 
  mutate(Species = gsub('C. optatus','Cuculus optatus', 
                        gsub('C. canorus','Cuculus canorus',Species_Latin)),
         LatJit = jitter(Latitude,amount = 1),
         LonJit = jitter(Longitude,amount = 1))

#and assign groups
all = all %>%
  mutate(Group = ifelse(Species == 'Cuculus canorus' & Latitude > 45 & Longitude < 5, 'CanWest',
                        ifelse(Species == 'Cuculus canorus' & Latitude > 35 & Longitude > 135 & Longitude < 165,'CanEast',
                               ifelse(Species == 'Cuculus optatus' & Latitude > 40 & Longitude > 75 & Longitude < 120,'OptWest',
                                      ifelse(Species == 'Cuculus optatus' & Longitude > 135 & Longitude < 150,'OptEast',
                                             'Other')))))

#and assign groups: VERY flexible, greater numbers 
all = all %>%
  mutate(Group = ifelse(Species == 'Cuculus canorus' & Latitude > 45 & Longitude < 40, 'CanWest',
                        ifelse(Species == 'Cuculus canorus' & Latitude > 35 & Longitude > 100 & Longitude < 165,'CanEast',
                               ifelse(Species == 'Cuculus optatus' & Longitude < 120,'OptWest',
                                      ifelse(Species == 'Cuculus optatus' & Longitude > 120 ,'OptEast',
                                             'Other')))))
all %>% count(Group)
alljit = all %>% mutate(loj = jitter(Longitude,amount=1),laj = jitter(Latitude,amount=1))
popgen_points = st_as_sf(alljit, coords = c("loj", "laj"), 
                         crs = 4326, agr = "constant", remove=F) 

#plot breeding polygons
cols = brewer.pal(12,'Paired')[c(1,2,5,6)]
map.polygon <- getMap(resolution = "low")
unrel_breed = ggplot() +
  geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='black',lwd=0.1) +
  geom_sf(data = popgen_points,aes(fill=Group,shape=Species),alpha=0.75,size=2,stroke=0.25)+
  #facet_wrap(Species~.)+
  xlab('') + ylab('') +
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=c(cols,'white'))+
  coord_sf(xlim = c(-10,185), 
           ylim = c(-30,75), expand = FALSE)+    #coord_equal() +  #FOR PUBLICATION MAIN FIGURE USE THIS, WON'T STRETCH / DISTORT LATITUDES
  theme_classic() +
  guides(fill=guide_legend(nrow=2,override.aes=list(shape=22)))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = 'top') +
  annotation_scale(line_width = 0.5)
unrel_breed

png('deprecated_sampling/double_sampling_2024AUG/Distribution_Groups_2024AUG12.png',units='in',res=300,height=6,width=10)
unrel_breed
dev.off()

all %>% filter(Group != 'Other') %>% group_by(Group) %>% sum_stats(MeanCoverage)

#retain n = 20 samples from each stronghold with the most similar coverage 
kept <- all %>% filter(Group != 'Other') %>% select(ID,Species,Group,Latitude,Longitude,MeanCoverage) %>% 
  mutate(Overall = mean(MeanCoverage), Diff = abs(MeanCoverage - Overall)) %>% # Calculate difference from mean  
  group_by(Species,Group) %>% slice_min(Diff,n=20)
kept %>% count(Species,Group)
covplot <- kept %>% group_by(Species,Group) %>% sum_stats(MeanCoverage) %>% 
  ggplot(aes(x=Group,y=mean,ymin=conf_low,ymax=conf_high,col=Group))+
  geom_errorbar()+
  geom_point(size=2.5)+ylab('Coverage by Group (mean & 95% CIs')+
  scale_color_manual(values=c(cols,'white'))+
  theme_bw()

png('deprecated_sampling/double_sampling_2024AUG/Coverage_Groups_2024AUG12.png',units='in',res=300,height=4,width=5)
covplot
dev.off()

kept %>% count(Group)

write.table(kept %>% ungroup %>% select(ID,Group),file='~/merondun/cuculus_migration/deprecated_sampling/double_sampling_2024AUG/Samples_Demography_N20_CCW-CCE-COW-COE_2024AUG12.pop',quote=F,sep='\t',row.names=F,col.names=F)
write.table(kept %>% ungroup %>% select(ID),file='~/merondun/cuculus_migration/deprecated_sampling/double_sampling_2024AUG/Samples_Demography_N20_CCW-CCE-COW-COE_2024AUG12.list',quote=F,sep='\t',row.names=F,col.names=F)


