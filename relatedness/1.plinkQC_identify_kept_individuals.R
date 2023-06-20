### Calculate relatedness
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/relatedness')
library(plinkQC)
library(tidyverse)

#first I ran plink, only on chromosome 6 for expediency, to examine relatedness
#~/modules/plink --geno 0.8 --mac 1 --allow-extra-chr --double-id --vcf ../merged/chr_6.vcf.gz --make-bed --out plink/chr_6

#run plink within R to calculate LD between SNPs. I will only look in 10 snp windows 
command <- "~/modules/plink --bfile plink/chr_6 --r2 --ld-window-kb 1000 --ld-window 10 --ld-window-r2 0.2 --out ld"
system(command)
df <- read.table('ld.ld', header=TRUE)

#sort the data by chromosome and position
df <- df[order(df$CHR_A, df$BP_A, df$BP_B),]

#save empty list to store the regions
regions <- list()

#start the first region with the first SNP pair
start <- df[1, 'BP_A']
end <- df[1, 'BP_B']
chr <- df[1, 'CHR_A']

#iterate over the SNP pairs
for (i in 2:nrow(df)) {
  #if the next SNP pair is on the same chromosome and overlaps with the current region, extend the region
  if (df[i, 'CHR_A'] == chr && df[i, 'BP_A'] <= end) {
    end <- max(end, df[i, 'BP_B'])
  } else {
    #ff the next SNP pair is on a different chromosome or does not overlap, save the current region and start a new one
    regions[[length(regions) + 1]] <- c(chr, start, end)
    start <- df[i, 'BP_A']
    end <- df[i, 'BP_B']
    chr <- df[i, 'CHR_A']
  }
}

#don't forget to save the last region
regions[[length(regions) + 1]] <- c(chr, start, end)

# Convert the regions into a data frame and save as a space-separated text file
regions_df <- do.call(rbind, regions)
colnames(regions_df) <- c('Chromosome', 'Start', 'End')
regions_df = regions_df %>% as.data.frame %>% mutate(ID = row_number())
write.table(regions_df, file='high_ld_regions.txt', row.names=FALSE, col.names=FALSE, sep=' ')

indir = 'plink'
qcdir = 'plink_qc'
name = 'chr_6'
path2plink = "/dss/dsshome1/lxc07/di39dux/modules/plink"
ldregions = 'high_ld_regions.txt'

fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
                                    path2plink=path2plink, 
                                    interactive=TRUE, verbose=TRUE,
                                    dont.check_sex=TRUE,
                                    do.run_check_ancestry = TRUE,
                                    do.evaluate_check_sex = FALSE,
                                    filter_high_ldregion=TRUE,
                                    high_ldregion_file=ldregions)

overview_individuals <- overviewPerIndividualQC(fail_individuals,
                                                interactive=TRUE)

Ids  <- cleanData(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
                  verbose=TRUE, showPlinkOutput=TRUE)

setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/allsamples/admixture')
library(tidyverse)
library(viridis)

#First, read in q matrices 
RUN='chr_10.AA-IF-GF_MM1-AC2-LD'
qs = list.files('.',paste0(RUN,'.*.Q$')) #find all files ending in .Q
famfile = paste0('../plink/',RUN,'.fam') #specify the .fam file 
samps = read.table(famfile,header=FALSE) %>% select(V2) %>% dplyr::rename(ID=V2)
qdat = NULL
for (q in qs){
  cat('Melting K: ',q,'\n')
  qf = read.table(q) #read in file
  k = ncol(qf) #k is simply the total number of columns
  names(qf) = paste0('K',seq(1,k)) #add 'KX' for all columns
  qfm = cbind(samps,qf) #bind with samples, ensure you're using the FAM file used for admixture! 
  qfm$MaxK = k   
  meltq = qfm %>% pivot_longer(!c(ID,MaxK),names_to = 'K',values_to = 'Q')  #melt it 
  qdat = rbind(qdat,meltq)
}

#read in metadata 
md = read.table('../EntireMetadata.txt',header=TRUE,comment.char = '',sep='\t')

#do a little PCA on distance
md1 = md %>%
  group_by(Species,Species_Latin) %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) %>%
  mutate(Migration = ifelse(grepl('AMBIG',Migration),'Ambiguous',
                            ifelse(grepl('north_south_one',Migration),'Ambiguous',
                                   ifelse(is.na(Migration),'unknown',Migration))),
         MigrationAlpha = ifelse(Migration == 'unknown',0.2,1))


#Prep
qdat = qdat %>% mutate(Kclust = paste0('K',MaxK))
qdat$Kclust <- factor(qdat$Kclust, levels = paste0("K", 2:10))
qm = left_join(qdat,md1)
ord = qm %>% select(ID,Distance) %>% unique %>% arrange(Distance)
qm$ID = factor(qm$ID,levels=ord$ID)

#Now plot
k2plot =
  qm %>% filter(MaxK < 7) %>%  #I only want to show K2-4
  ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
  geom_col(color = NA) +
  facet_grid(Kclust ~ Species_Latin, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "",y = "Ancestry Coefficient") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_viridis('K',discrete=TRUE,option='turbo')+ 
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(size=5,angle=90,vjust=0,hjust=0),
    panel.grid = element_blank(),
    strip.text = element_text(size = 5),
    legend.position='top')
k2plot

pdf('Admixture_AllSamples_n408.pdf',height=4,width=9)
png('Admixture_AllSamples_n408.png',height=4,width=9,units='in',res=600)
k2plot
dev.off()



