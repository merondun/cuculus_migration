### Plot ADMIXTURE
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/full_qs/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(terra)
library(tidyverse)
library(viridis)
library(stringr)
library(meRo) #devtools::install_github('merondun/meRo')
library(LEA) #install with bioconductor, you don't actually need this if you impot your own q-matrix
library(tess3r) #install with github devtools
library(rworldmap) #for ggplot mapping 
library(sf) #for spatial plotting of distance vectors to confirm
library(ggspatial) #to add scale bars onto maps
library(ggpubr)

prefixes = gsub('.2.Q','',list.files('.',pattern='.*.2.Q'))
prefixes = prefixes[grepl('^Canorus.MQ|^Optatus.MQ',prefixes)]
prefixes
qdir = '.' #directory with Q files
counter = 0 

for (admix_run in prefixes) { 
  counter = counter + 1
  cat('Working on run: ',admix_run,'\n')
  prefix = admix_run 
  admix = melt_admixture(prefix = prefix, qdir = qdir)
  
  #read in metadata
  md = read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')
  admixmd = left_join(admix,md)
  
  #Count number of SNPs used 
  nsnps = str_extract(readLines(paste0(qdir,'/',admix_run,'.log2.out')), "(?<=Size of G: )\\S+") %>% na.omit
  nsnps = nsnps[[1]]
  snplab = sprintf('n = %s, n = %.2fM SNPs', unlist(strsplit(nsnps, "x"))[1], as.numeric(unlist(strsplit(nsnps, "x"))[2]) / 10^6)
  
  #Make a readable label 
  label = gsub('\\..*','',admix_run)
  cat(label,' for run: ',admix_run,', plot: ',paste0('p',counter),'\n')
  
  #add a distance vector based on lat/long PC1
  admixmd = admixmd %>% 
    mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
    { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
    rename(Distance = pca) 
  
  admixmd = admixmd %>% mutate(ID = fct_reorder(ID,desc(Distance)),
                               Species = gsub('CB','CC',Species),
                               Species_Latin = gsub('C. canorus bakeri','C. canorus',Species_Latin),
                               Species_Latin = gsub('C. canorus canorus','C. canorus',Species_Latin)) 
  admixmd %>% filter(ID != '445_CC_UNK_TOG_F') #remove the individual museum sample from Africa for tesselation, which we don't have 
  #line number 231 from the fam file 
  
  #loop through all admixture runs and extract the average correlation values from evalAdmix, we want to MINIMIZE this! (closest to 0)
  evaldat = NULL; for (Kval in seq(2,10,1)){
    r <- as.matrix(read.table(paste0("eval_",admix_run,'_',Kval)))
    mean_value <- mean(r,na.rm=TRUE)
    median_value <- median(r,na.rm=TRUE)
    sd_value <- sd(r,na.rm=TRUE)
    iqr_value <- IQR(r,na.rm=TRUE)
    valdat = data.frame(K = Kval,mean = mean_value,median=median_value,sd=sd_value,iqr=iqr_value)
    evaldat = rbind(valdat,evaldat)
  }
  
  #plot, for main figure show the n=3 lowest median
  targs = evaldat %>% slice_min(abs(median),n=3)
  ep = evaldat %>% 
    ggplot(aes(x=K,y=median,ymin=median-iqr,ymax=median+iqr))+
    geom_rect(data=targs,aes(xmin=K-0.25,xmax=K+0.25,ymin=-Inf,ymax=Inf),fill='darkseagreen3')+
    geom_text(aes(y = 0.014, label = format(signif(median, 2), scientific = TRUE)),size=2) +  ylim(c(-0.015,0.015))+
    geom_hline(yintercept=0,lty=2)+
    geom_point()+ylab('Median +/- IQR Correlation of Residuals') +
    geom_errorbar()+
    theme_bw() + 
    ggtitle(paste0(label,': ',snplab))+
    scale_x_continuous(breaks = seq(min(evaldat$K), max(evaldat$K), by = 1)) +
    coord_flip()
  ep
  assign(paste0('e',1),ep)
  
  png(paste0('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/evalAdmix_',admix_run,'_AllK_2024MAR04.png'),res=300,units='in',height=5,width=6)
  print(e1)
  dev.off()
  
  #Now plot ADMIXTURE 
  adplot =
    admixmd %>% filter(Specified_K <= 5 & Specified_K != 4) %>%  #specify the levels you want 
    mutate(Specified_K = paste0('K',Specified_K)) %>% 
    ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
    geom_col(color = "gray", size = 0.1) +
    facet_grid(Specified_K~Species, scales = "free", space = "free") +
    theme_minimal(base_size=6) + labs(x = "",y = "") +
    scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
    scale_fill_viridis(discrete=TRUE,option='turbo')+
    ggtitle(paste0(label,'\n', snplab))+
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=3),
      panel.grid = element_blank(),
      legend.position = 'bottom',
      plot.title = element_text(size=6)
    )
  adplot
  assign(paste0('p',counter),adplot)
  
  ### Tesselation
  #Specify K for plotting on map
  show_k = 2
  show_q = read.table(paste0(qdir,'/',admix_run,'.',show_k,'.Q')) #read in the specific file
  show_q_mat = as.matrix(show_q) #convert it to a matrix
  class(show_q_mat) = c('tess3Q','matrix','array') #make sure tess3r thinks that it's actually a tess object
  coords = admixmd %>% select(ID,Longitude,Latitude) %>% unique %>% select(Longitude,Latitude) #convert lat and long
  coords_mat = as.matrix(coords) #convert coordinates to matrix

  #plot using their base function
  plot(show_q_mat, coords_mat, method = "map.max", interpol = FieldsKrigModel(10),
       main = "Ancestry coefficients",
       xlab = "Longitude", ylab = "Latitude",
       resolution = c(300,300), cex = .4,
       col.palette = viridis(show_k))

  #plot using ggplot
  map.polygon <- getMap(resolution = "low")
  pl = ggtess3Q(show_q_mat, coords_mat, map.polygon = map.polygon,col.palette = viridis(show_k))
  pl +
    geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
    xlim(min(coords$Longitude)-5,max(coords$Longitude)+5) +
    ylim(min(coords$Latitude)-5,max(coords$Latitude)+5) +
    coord_equal() +
    geom_point(data = coords, aes(x = Longitude, y = Latitude), size = 2,pch=21,fill='white') +
    xlab("Longitude") + ylab("Latitude") + ggtitle(paste0('K = ',show_k))+
    theme_classic()

  #plot all K...
  maximum_k = 5
  map.polygon <- getMap(resolution = "low")
  for (kval in seq(2,maximum_k,1)) {

    cat('Working on K = ',kval,'\n')
    show_q = read.table(paste0(qdir,'/',admix_run,'.',kval,'.Q')) #read in the specific file
    if (label != 'Optatus') {
      fam = read.table(paste0(admix_run,'.fam'))
      af_ind = which(fam$V2 == "445_CC_UNK_TOG_F")
      show_q = show_q[-af_ind,] #for canorus, exclude the african individual
    } else { cat('Not removing African individual, optatus\n')}

    show_q_mat = as.matrix(show_q) #convert it to a matrix
    class(show_q_mat) = c('tess3Q','matrix','array') #make sure tess3r thinks that it's actually a tess object
    coords = admixmd %>% select(ID,Longitude,Latitude) %>% filter(!grepl('445_CC_UNK_TOG_F',ID)) %>% unique %>% select(Longitude,Latitude) #convert lat and long
    coords_mat = as.matrix(coords) #convert coordinates to matrix

    #plot
    pl = ggtess3Q(show_q_mat, coords_mat, map.polygon = map.polygon,col.palette = viridis(kval))
    pp = pl +
      geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
      xlim(min(coords$Longitude)-5,max(coords$Longitude)+5) +
      ylim(min(coords$Latitude)-5,max(coords$Latitude)+5) +
      #coord_equal() +  #FOR PUBLICATION MAIN FIGURE USE THIS, WON'T STRETCH / DISTORT LATITUDES
      geom_point(data = coords, aes(x = Longitude, y = Latitude), size = 2,pch=21,fill='white') +
      xlab("Longitude") + ylab("Latitude") + ggtitle(paste0('K = ',kval))+
      theme_classic()
    
    kept = admixmd %>% select(ID,Specified_K,K,Q,Latitude,Longitude,Group = Cluster) %>% 
      filter(!grepl('445_CC_UNK_TOG_F',ID))
    group_summaries = kept %>% 
      filter(Specified_K == kval) %>%
      group_by(Group,K) %>%  #within each group and K, average lat/long/q and count number of individuals 
      summarize(Lat = mean(Latitude),
                Long = mean(Longitude),
                Q = mean(Q),
                N = n_distinct(ID)) %>% 
      ungroup() %>% #
      #Calculate scaling factors for the pies based on num samples
      mutate(MinN = min(N),
             MaxN = max(N)) %>%
      group_by(Group) %>%
      mutate(Scaling_factor = ((N - MinN) / (MaxN - MinN) * 10) + 2) %>%
      select(-MinN, -MaxN) 
    
    ##### Plot pies across the world 
    plot_pie <- function(data) {
      ggplot(data, aes(x = "", y = Q, fill = K,)) +
        geom_bar(col='white',lwd=0.5,width = 1, stat = "identity") +
        coord_polar("y") +
        scale_fill_viridis(discrete=TRUE)+
        theme_void() +
        theme(legend.position = "none")
    }
    
    #set up map and make a sf object from the summaries 
    sites = st_as_sf(group_summaries, coords = c("Long", "Lat"), crs = 4326, agr = "constant") 
    
    # Main map plot
    p = 
      # ggplot()+ #for showing labels only 
      # geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
      ggtess3Q(show_q_mat, coords_mat, map.polygon = map.polygon,col.palette = viridis(kval)) + 
      geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
      geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1) +
      xlab('')+ylab('')+
      coord_sf(xlim = c(min(group_summaries$Long)-5, max(group_summaries$Long)+5), 
               ylim = c(min(group_summaries$Lat)-5, max(group_summaries$Lat)+5), expand = FALSE)+
      theme_classic(base_size = 8)+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
      theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
    p
    
    # Add pies
    for (i in unique(group_summaries$Group)) {
      subset_data = group_summaries %>% filter(Group == i)
      lon = unique(subset_data$Long)
      lat = unique(subset_data$Lat)
      scale_factor = unique(subset_data$Scaling_factor)
      cat('Scaling factor is : ',scale_factor,' for group : ',i,'\n')
      pie = plot_pie(subset_data)
      p <- p + annotation_custom(ggplotGrob(pie), 
                                 xmin = lon - scale_factor, 
                                 xmax = lon + scale_factor, 
                                 ymin = lat - scale_factor, 
                                 ymax = lat + scale_factor)
    }
    p
    assign(paste0('t',kval),p)

  }

  png(paste0('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Tesselation_',admix_run,'_AllK_2024MAR04.png'),res=300,units='in',height=7,width=10)
  print(ggarrange(t2,t3,t4,t5,ncol=2,nrow=2))
  dev.off()
  
}

png('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Admixture_Sensitivity_Migration-CanOpt_2024MAR06.png',height=4,width=8,units='in',res=600)
ggarrange(p1,p4,nrow=2,ncol=1,common.legend = TRUE)
dev.off()

##### This next section simply adds 'K1/K2' values for K=2 to assign eastern/western groups. 
#Grab metadata sheet for canorus and optatus separately
mdall = NULL # for the final metadata sheet with canorus and optatus separately

#read in metadata
md = read.table('~/merondun/cuculus_migration/Full_Metadata.txt',header=TRUE,comment.char='',sep='\t') %>% as_tibble

#read in can and opt 
can = "Autosomes.Canorus.IF-GF-MM1-AA-BP-NEUTRAL-P25.10.R1-UNLINKED"
canmix = melt_admixture(prefix = can, famdir = famdir, qdir = qdir)
canmixmd = left_join(canmix,md)

#add a distance vector based on lat/long PC1
canmixmd = canmixmd %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#to simply write a wide-format metadata for canorus and optatus 
cansub = canmixmd %>% filter(Specified_K == '2') %>% pivot_wider(names_from = K, values_from = Q)
cansub = cansub %>% mutate(Population = ifelse(K1 < 0.2,'Canorus_East',
                                               ifelse(K1 > 0.8,'Canorus_West','Intermediate')))
cansub %>% ggplot(aes(x=Longitude,y=Latitude,col=Population))+geom_point()+theme_bw()

#read in opt 
opt = "Autosomes.Optatus.IF-GF-MM1-AA-BP-NEUTRAL-P25.10.R1-UNLINKED"
optmix = melt_admixture(prefix = opt, famdir = famdir, qdir = qdir)
optmixmd = left_join(optmix,md)

#add a distance vector based on lat/long PC1
optmixmd = optmixmd %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 

#to simply write a wide-format metadata for optorus and optatus 
optsub = optmixmd %>% filter(Specified_K == '2') %>% pivot_wider(names_from = K, values_from = Q)
optsub = optsub %>% mutate(Population = ifelse(K1 < 0.2,'Optatus_West',
                                               ifelse(K1 > 0.8,'Optatus_East','Intermediate')))
optsub %>% ggplot(aes(x=Longitude,y=Latitude,col=Population))+geom_point()+theme_bw()

#bind together
mdall = rbind(cansub,optsub)
write_tsv(mdall,file='../Populations_CanorusOptatus_K2__2023AUG25.txt')


#plot pie charts 
admix = melt_admixture(prefix = 'Optatus.MQ-5X-MM1-AA-LDr2w50', qdir = qdir)

#read in metadata
md = read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')
admixmd = left_join(admix,md)
kept = admixmd %>% select(ID,Specified_K,K,Q,Latitude,Longitude,Group = Cluster)
kept
#columns ID, specified K (e.g. the maximum K), K (inferred cluster within each K level), Q, lat/long, group 
# # A tibble: 12,474 × 7
# ID               Specified_K K           Q Latitude Longitude Group
# <chr>                  <dbl> <chr>   <dbl>    <dbl>     <dbl> <chr>
#   1 001_CB_ORW_CHN_F          10 K1    0.00001     47.5      125. CC_1 
# 2 001_CB_ORW_CHN_F          10 K2    0.00001     47.5      125. CC_1 
# 3 001_CB_ORW_CHN_F          10 K3    0.00001     47.5      125. CC_1 
# 4 001_CB_ORW_CHN_F          10 K4    0.00001     47.5      125. CC_1

#count the proportion of each haplogroup within each distance group 
plot_k = 3 #change this to the K value you want to plot 
group_summaries = kept %>% 
  filter(Specified_K == plot_k) %>%
  group_by(Group,K) %>%  #within each group and K, average lat/long/q and count number of individuals 
  summarize(Lat = mean(Latitude),
            Long = mean(Longitude),
            Q = mean(Q),
            N = n_distinct(ID)) %>% 
  ungroup() %>% #
  #Calculate scaling factors for the pies based on num samples
  mutate(MinN = min(N),
         MaxN = max(N)) %>%
  group_by(Group) %>%
  mutate(Scaling_factor = ((N - MinN) / (MaxN - MinN) * 10) + 2) %>%
  select(-MinN, -MaxN) 

##### Plot pies across the world 
plot_pie <- function(data) {
  ggplot(data, aes(x = "", y = Q, fill = K,)) +
    geom_bar(col='black',lwd=0.5,width = 1, stat = "identity") +
    coord_polar("y") +
    scale_fill_viridis(discrete=TRUE)+
    theme_void() +
    theme(legend.position = "none")
}

#set up map and make a sf object from the summaries 
world = map_data("world")
sites = st_as_sf(group_summaries, coords = c("Long", "Lat"), crs = 4326, agr = "constant") 

# Main map plot
p = 
  ggplot()+ #for showing labels only 
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(group_summaries$Long)-5, max(group_summaries$Long)+5), 
           ylim = c(min(group_summaries$Lat)-5, max(group_summaries$Lat)+5), expand = FALSE)+
  theme_classic(base_size = 8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
p

# Add pies
for (i in unique(group_summaries$Group)) {
  subset_data = group_summaries %>% filter(Group == i)
  lon = unique(subset_data$Long)
  lat = unique(subset_data$Lat)
  scale_factor = unique(subset_data$Scaling_factor)
  cat('Scaling factor is : ',scale_factor,' for group : ',i,'\n')
  pie = plot_pie(subset_data)
  p <- p + annotation_custom(ggplotGrob(pie), 
                             xmin = lon - scale_factor, 
                             xmax = lon + scale_factor, 
                             ymin = lat - scale_factor, 
                             ymax = lat + scale_factor)
}
p

png('~/merondun/cuculus_migration/admixture/Optatus_K3_SpatialQ_2024MAR06.png',units='in',res=600,height=3.5,width=7)
p
dev.off()

