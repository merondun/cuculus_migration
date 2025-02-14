#### Determine relatedness with PHI statistic 
setwd('~/merondun/cuculus_migration/relatedness/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
#Igraph approach
library(tidyverse)
library(RColorBrewer)
library(geosphere)
library(igraph)
library(spThin)
library(sf)
library(ggspatial)
library(factoextra)
library(ggpubr)
library(viridis)

# Read metadata and filter for necessary columns
md = read_tsv('../Full_Metadata.txt') 
species = c('CC','CO')
world <- map_data("world")
set.seed(9999)

#Find samples within specified distance 
assign_clusters_based_on_distance <- function(ds) {
  #Calculate pairwise distances
  distances <- geosphere::distm(ds[, c("Longitude", "Latitude")])
  
  # Create an adjacency matrix where 1 represents distances within 500km
  adjacency_matrix <- distances <= 500000
  
  # Create a graph frm the adjacency matrix
  g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", diag = FALSE)
  
  # Find connected clusters
  clusters <- components(g)$membership
  
  return(clusters)
}

#apply to each species 
ds_clustered <- md %>% 
  filter(grepl('CC|CO',SpeciesShort)) %>% 
  group_by(SpeciesShort) %>%
  mutate(cluster = assign_clusters_based_on_distance(cur_data())) %>%
  ungroup()


dist_jit = ds_clustered %>% mutate(Cluster = paste0(SpeciesShort,'_',cluster)) 
colshape = dist_jit %>% ungroup %>% select(Cluster) %>% unique %>% 
  mutate(ClustCol = rep_len(c(viridis(12, option = 'turbo'),brewer.pal(12, 'Paired')), 35), ClustShape = rep_len(c(21, 24, 25, 22,23), 35)) 

dist_site = st_as_sf(dist_jit, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") 
distp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey50',fill='grey95') +
  geom_sf(data = dist_site, 
          aes(fill=Cluster,shape=Cluster),
          size=3,show.legend = T) +
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_manual(values=colshape$ClustCol,breaks=colshape$Cluster)+
  scale_shape_manual(values=colshape$ClustShape,breaks=colshape$Cluster)+
  coord_sf(xlim = c(min(ds_clustered$Longitude)-5, max(ds_clustered$Longitude)+5), 
           ylim = c(min(ds_clustered$Latitude)-5, max(ds_clustered$Latitude)+5), expand = FALSE)+
  theme_classic()+
  facet_grid(SpeciesShort ~ .)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)
distp

pdf('~/merondun/cuculus_migration/relatedness/DistanceClusters_2024FEB28.pdf',height=6,width=7)
distp
dev.off()

write.table(dist_jit %>% ungroup %>% select(ID,Cluster),'~/merondun/cuculus_migration/relatedness/DistanceClusters_2024FEB28.txt',quote=F,sep='\t',row.names=F)
