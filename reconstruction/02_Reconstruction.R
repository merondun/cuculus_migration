#### Plot & estimate species tree 
setwd('~/EvoBioWolf/CUCKOO_migration/reconstruction/')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ggtree)
library(ape)
library(tidyverse)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(phytools)
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(magick)

# # Only run this the first time to convert the nex trees
# # Read in nexus
# trees <- read.nexus('hackett_sequenced_species.nex')
# 
# # Re-write trees
# write.tree(trees, file = "hackett_sequenced_species.tre")

# Read in iqtree
iqtree <- read.iqtree('species.tre')
iqtree_root <- treeio::root(as.phylo(iqtree),outgroup='Geococcyx_californianus',resolve.root=TRUE)

gg <- ggtree(iqtree_root, layout = "rectangular")
gg
gg$data <- gg$data %>% mutate(label = gsub('.*\\/','',label))
gg + 
  geom_tippoint()+
  geom_tiplab(align = TRUE)+
  #xlim(c(min(gg$data$x)-.05,max(gg$data$x)*1.3))+
  geom_nodelab(geom='text',mapping=aes(label=label),col='black',size=3,show.legend=F,vjust=-1,hjust=1.15)+
  geom_treescale(y=10,x=0.01)+
  theme(legend.position = 'none')


# Generate base tree 
md <- read.table('Species_Data.txt',header=TRUE)
targ_tree <- as.phylo(iqtree_root)
ggt <- ggtree(targ_tree, layout = "rectangular",branch.length='none') %<+% md 

#grab only residency
phenos <- as.data.frame(ggt$data %>% filter(isTip == TRUE))
mat <- as.matrix(phenos %>% select(Residency))
phenotypes <- setNames(as.vector(mat[, 1]), targ_tree$tip.label)

#Plot probabilities 
t2 <- multi2di(targ_tree)
t2$edge.length[is.na(t2$edge.length)] <- 1e-8

# Quick ape method 
fitER <- ape::ace(phenotypes,t2,model="ER",type="discrete")

# Extract nodes and the proportions for pies
nodes <- data.frame(
  node=1:t2$Nnode+Ntip(t2),
  fitER$lik.anc)

t2$node.label <- gsub('.*\\/','',t2$node.label)
tip_rename <- left_join(data.frame(label=t2$tip.label),md) %>% 
  mutate(lab = paste0(Genus,'_',Species))
t2$tip.label <- tip_rename$lab

#### Main plot, node labels < 100 bootstrap support 
pdf('20250115_ReconstructionER-Cuculiformes_Residency.pdf',height=7,width=7)
cols<-setNames(viridis(3)[1:length(unique(phenotypes))],sort(unique(phenotypes)))
plotTree(t2,ftype="off")

#Node points
node_labs <- ifelse(round(as.numeric(t2$node.label)) == 100, NA, round(as.numeric(t2$node.label)))
nodelabels(node_labs,node=1:t2$Nnode+Ntip(t2),cex=0.7,bg=NA,adj = c(-1))
nodelabels(node=1:t2$Nnode+Ntip(t2),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.5,
                  y=max(nodeHeights(t2)-50),fsize=2)
dev.off()

##### Supp plot, show tip points, no node labels 
png('20250115_ReconstructionER-Cuculiformes_Residency-TIPS.png',height=7,width=7,units='in',res=300)
cols<-setNames(viridis(3)[1:length(unique(phenotypes))],sort(unique(phenotypes)))
plotTree(t2,ftype="off")

# Tip points 
tip_colors <- cols[phenotypes] 
tiplabels(pch=24, bg=tip_colors, col="black", cex=1, adj = c(2.5, 0.5))

#Node points
node_labs <- ifelse(round(as.numeric(t2$node.label)) == 100, NA, round(as.numeric(t2$node.label)))
nodelabels(node=1:t2$Nnode+Ntip(t2),
           pie=fitER$lik.anc,piecol=cols,cex=0.4)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.5,
                  y=max(nodeHeights(t2)-50),fsize=2)
dev.off()


#### Tip labels 
pdf('20250115_ReconstructionER-Cuculiformes_Residency-TIPLABELS.pdf',height=7,width=9)
cols<-setNames(viridis(3)[1:length(unique(phenotypes))],sort(unique(phenotypes)))
plotTree(t2,fsize=0.8,ftype="i")
tiplabels(pie=to.matrix(phenotypes,sort(unique(phenotypes))),piecol=cols,cex=0.05)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.5,
                  y=max(nodeHeights(t2)-50),fsize=2)
dev.off()

##### Migration route reconstruction ####
#grab only route
mat2 <- as.matrix(phenos %>% select(Route))
phenotypes2 <- setNames(as.vector(mat2[, 1]), targ_tree$tip.label)

# Quick ape method 
fitER2 <- ape::ace(phenotypes2,t2,model="ER",type="discrete")

# Extract nodes and the proportions for pies
nodes2 <- data.frame(
  node=1:t2$Nnode+Ntip(t2),
  fitER2$lik.anc)

#### Plot nodes and node labels 
pdf('20250115_ReconstructionER-Cuculiformes_Route.pdf',height=7,width=7)
cols<-setNames(viridis(7,option='turbo')[1:length(unique(phenotypes2))],sort(unique(phenotypes2)))
plotTree(t2,ftype="off")

node_labs <- ifelse(round(as.numeric(t2$node.label)) == 100, NA, round(as.numeric(t2$node.label)))
nodelabels(node_labs,node=1:t2$Nnode+Ntip(t2),cex=0.7,bg=NA,adj = c(-1))
nodelabels(node=1:t2$Nnode+Ntip(t2),
           pie=fitER2$lik.anc,piecol=cols,cex=0.4)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=max(nodeHeights(t2)-50),fsize=2)
dev.off()

#### Plot TIP SHAPES 
png('20250115_ReconstructionER-Cuculiformes_Route-TIPS.png',units='in',res=300,height=7,width=7)
cols<-setNames(viridis(7,option='turbo')[1:length(unique(phenotypes2))],sort(unique(phenotypes2)))
plotTree(t2,ftype="off")

# Tip points 
tip_colors2 <- cols[phenotypes2] 
tiplabels(pch=24, bg=tip_colors2, col="black", cex=1, adj = c(2.5, 0.5))

node_labs <- ifelse(round(as.numeric(t2$node.label)) == 100, NA, round(as.numeric(t2$node.label)))
nodelabels(node=1:t2$Nnode+Ntip(t2),
           pie=fitER2$lik.anc,piecol=cols,cex=0.4)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=max(nodeHeights(t2)-50),fsize=2)
dev.off()


# Tip labels 
pdf('20250115_ReconstructionER-Cuculiformes_Route-TIPLABELS.pdf',height=7,width=7)
cols<-setNames(viridis(7,option='turbo')[1:length(unique(phenotypes2))],sort(unique(phenotypes2)))
plotTree(t2,fsize=0.8,ftype="i")
tiplabels(pie=to.matrix(phenotypes2,sort(unique(phenotypes2))),piecol=cols,cex=0.1)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=max(nodeHeights(t2)-50),fsize=2)
dev.off()


# For plotting with ggtree 
## cols parameter indicate which columns store stats
residency_colors <- data.frame(cols=viridis(3),levs=names(nodes[,2:4]))
pies <- nodepie(nodes, cols=2:4,outline.color='black',outline.size = 0.1)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = residency_colors$cols,breaks=residency_colors$levs))

t3 <- full_join(t2, data.frame(label = names(phenotypes), stat = phenotypes ), by = 'label')
tp <- ggtree(t3,layout='rectangular') %<+% md
tp$data$dummy <- 1
tp_final <- tp + geom_inset(pies, width = .09, height = .09)
tp_phenos <- tp_final +
  geom_tippoint(aes(fill=Residency),pch=21,size=1.5)+
  scale_fill_manual(values = residency_colors$cols,breaks=residency_colors$levs)
tp_phenos

ggsave('20250115_Cuculiformes-Residency-ggtree.pdf',tp_phenos,height=15,width=15,dpi=300)

