setwd('F:/Research/scratch/cuckoo/2022_04/QC/alignments/')
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggsci)
library(ggpubr)
library(viridis)
library(openxlsx)

#import alignments
aln <- read.table('Alignments.txt',header=T)
aln <- aln %>% rename(ID = SAMPLE)
coord <- read.xlsx('D:/Sync/JM/Research/Cuckoo/Cuculus_Metadata_n404.xlsx')
alnm  <- merge(aln,coord,by='ID')
raw <- read.table('reads.txt',header=T)
raw <- raw %>% group_by(ID) %>% summarize(total=sum(Reads))
ra <- merge(alnm,raw,by='ID')
ra$aligned <- ra$TOTAL_READS/1000000
ra$rate <- ra$aligned/ra$total
cs1 <- subset(ra,Species != 'EC' & Sex != 'U' & Species != 'HH')
cs1 <- cs1 %>% mutate(Species = gsub('CC','C. canorus canorus',Species),
                      Species = gsub('CB','C. canorus bakeri',Species),
                      Species = gsub('CO','C. optatus',Species),
                      Species = gsub('CP','C. poliocephalus',Species),
                      Species = gsub('CM','C. micropterus',Species),
                      Species = gsub('CS','C. saturatus',Species))
cs1$Species <- factor(cs1$Species,levels=c('C. canorus canorus','C. canorus bakeri','C. optatus','C. saturatus','C. micropterus','C. poliocephalus'))
ccol <- coord %>% select(DistanceClade,DistanceColor) %>% unique()
cp <- cs1 %>% ggplot(aes(x=Structure_Order,y=rate,col=DistanceClade,shape=Species))+
  geom_point()+
  scale_color_manual(values=col$DistanceColor,
                     breaks=col$DistanceClade)+
  scale_shape_manual(values=c('C. canorus canorus'=16,'C. canorus bakeri'=15,'C. optatus'=17,'C. saturatus'=6,'C. micropterus'=3,'C. poliocephalus'=4))+
  facet_grid(.~Sex,scale='free')+
  ylab('Mean Coverage')+xlab('Sample (Same Order as Structure Plots)')+
  theme_minimal()+
  labs(color='K-Means Distance Clade') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
cp
cs1 %>% group_by(Species) %>% summarize(mean=mean(rate))
png('Alignments_Cuckoo_All.png',height=7,width=10,res=300,bg='white',units='in')
cp
dev.off()

