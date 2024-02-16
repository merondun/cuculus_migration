### Calculate relatedness
.libPaths('~/mambaforge/envs/r/lib/R/library')  # Setting the path to the library
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/relatedness')  # Setting the working directory
library(tidyverse)  # Loading tidyverse library for data manipulation

# Reading relatedness data and renaming columns
r = read.table('chr_10.rel.input',header=TRUE) %>% select(INDV1,INDV2,RELATEDNESS_PHI) %>% as_tibble
names(r) = c('SampleA','SampleB','phi')

# Reading missingness data, renaming columns, and merging with relatedness data
m = read.table('chr_10.miss.input',header=TRUE) %>% select(INDV,F_MISS) %>% as_tibble
names(m) = c('SampleA','missA')
rm = left_join(r,m) %>% 
  left_join(.,m %>% dplyr::rename(SampleB=SampleA,missB=missA))

# Filtering samples with a certain relatedness threshold, done iteratively
choice = rm %>% filter(SampleA != SampleB & phi > 0.0884)
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
#total kept
samples_to_keep <- samples[!(samples %in% removed[, 1])]
write.table(samples_to_keep,file='Unrelated_Samples_2023JUNE21.txt',quote=F,sep='\t',row.names=F,col.names=F)  # Writing kept samples to a file

