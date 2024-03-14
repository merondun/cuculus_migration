### Remove relatives
.libPaths('~/mambaforge/envs/r/lib/R/library')  # Setting the path to the library
setwd('~/merondun/cuculus_migration/relatedness/')  # Setting the working directory
library(tidyverse)  # Loading tidyverse library for data manipulation

# Reading relatedness data and renaming columns
r = read_tsv('PLINK_Relatedness_2024FEB29.txt') 
names(r) = c('SampleA','SampleB','phi')
r = r %>% filter(SampleA != SampleB) %>% unique

# Reading missingness data, renaming columns, and merging with relatedness data
m = read_tsv('missingness.txt') %>% unique %>% select(INDV,N_GENOTYPES_FILTERED,F_MISS)
names(m) = c('Sample','Sites','Missing')
phi = left_join(r,m %>% select(SampleA = Sample, missA = Missing)) %>% 
  left_join(.,m %>% select(SampleB = Sample, missB = Missing))

#add more phenotypes
md = read_tsv('../Full_Metadata.txt')
phi = left_join(phi,md %>% select(SampleA = ID, migA = Migration)) %>% 
  left_join(.,md %>% select(SampleB = ID, migB = Migration))

# Filtering samples with a certain relatedness threshold, done iteratively
choice = phi
samples = unique(c(choice$SampleA,choice$SampleB)) ; length(unique(samples))
removed = NULL
for (run in seq(1,5,1)) { 
  
  cat('Running ',run,'\n')
  # Iterating over each sample
  for (samp in samples) {
    # Grabbing all the records with the current sample
    sf <- choice %>% filter(samp == SampleA | samp == SampleB) %>% arrange(desc(phi))
    if(nrow(sf) < 1) next
    # Deciding which individual to remove based on missing data
    bad <- sf %>% mutate(Remove = ifelse(missA > missB, SampleA,
                                         ifelse(missB > missA, SampleB,
                                                SampleA))) %>%
      pull(Remove) %>% head(n=1)
    
    # Unless there are no samples to remove, continue to the next iteration
    if(length(bad) < 1) next
    # Removing the bad sample from the pool
    cat('Removing bad sample: ',bad,'\n')
    choice = choice %>% filter(SampleA != bad & SampleB != bad)
    # And then restarting the whole process
    removed <- rbind(removed,bad)
    rm(bad)
  }
  
}

#total removed
length(unique(removed))
write.table(removed,file='Removed_Samples_2024FEB29.txt',quote=F,sep='\t',row.names=F,col.names=F)  # Writing related samples to file 


