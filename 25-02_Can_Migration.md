#  Cuculus Speciation: Migratory Biogeographic Genomics

```bash
title: "Cuculus Migratory Genomics" 
author: "Justin Merondun"
date: "`r Sys.Date()`"
output:
  html_document:
    theme: journal
    toc: true
    toc_float:
      collapsed: true
```

## Summarize Metadata

```R
setwd('C:/Users/herit/Dropbox/CUCKOO_migration/Manuscript/Supplement/SupplementaryTables/')
library(tidyverse)
library(readxl)
library(meRo)

# Load the Excel file, skipping the first row
df <- read_excel("Supplementary_Tables.xlsx", skip = 1,sheet = 'S1') %>% as_tibble()

df %>% count(SpeciesShort)
df %>% mutate(Source = ifelse(Source == 'Blood','Blood',
                              ifelse(grepl('Feather',Source),'Feather',
                                     'Muscle'))) %>% 
  count(Source)
summary(df$Raw_Bases_Gb)
sum(df$Raw_Bases_Gb)
sum(df$Raw_Reads_M)
summary(df$Mapped_Reads_M)
summary(df$MeanCoverage)
```



## Trimming

Business as usual; align and trim SRRs. 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00
RUN=$1
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/rawdata/WGS/ILLUMINA/
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq

RUN=$1

file1=$(ls ${wd}/${RUN}__R1__*)
file2=$(ls ${wd}/${RUN}__R2__*)

SCRATCH=/tmp/$SLURM_JOB_ID

#adapter trim
bbduk.sh t=6 -Xmx24g in=${file1} in2=${file2} out=${outdir}/${RUN}.trim.fastq.gz minlen=25 qtrim=rl trimq=2 ktrim=r k=23 mink=11 ref=~/modules/bbmap/adapters.fa hdist=1 tpe tbo
```

and submit:

```bash
for i in $(cat Libraries.list); do sbatch -J $i Trim.sh $i; done 
```

## Identify Read Group

This script will find read groups from the fastq files, provided you have normal fastq, and your files are named in our convention (SM__SRR.something.fastq.gz).

```bash
#!/bin/bash

#In a directory with the cleaned .fq.gz files, this will create ID (for us just SRR ID), SM (most important, for VCF files), LB (in case same library across multiple lanes), and PU (run data) fields. It will also count the number of SRRs per SM for merging later. 

mkdir RGs

for i in $( ls *.trim.fastq.gz | sed 's/\..*//g' ); do
        zcat $i.trim.fastq.gz | head -n 1 > ./RGs/$i.txt
        awk 'BEGIN {FS=":"}{print FILENAME, $3, $4, $NF}' ./RGs/$i.txt > ./RGs/$i.tmp
        ex -sc '%s/.txt//g' -c x ./RGs/$i.tmp
        ex -sc '%s:./RGs/::g' -c x ./RGs/$i.tmp

        #okay, now subet the parts we want
        sed 's/__.*//g' ./RGs/$i.tmp > ./RGs/$i.SM
        sed 's/.*__//g' ./RGs/$i.tmp | cut -d ' ' -f 1 > ./RGs/$i.ID
        awk '{print $4}' ./RGs/$i.tmp > ./RGs/$i.tmp2
        paste ./RGs/$i.SM ./RGs/$i.tmp2 | sed 's/\t/./g' > ./RGs/$i.LB

        #run details for PU
        zcat $i.trim.fastq.gz |head -n 1 |cut -f 3-5 -d ":" > ./RGs/$i.PU

        #Determine number of runs per sample, and output to .numlibs file
        ID="$(cat ./RGs/${i}.ID | tail -1)";
        SM="$(cat ./RGs/${i}.SM | tail -1)";
        echo ${ID} >> ./RGs/${SM}.libco
        cat ./RGs/${SM}.libco | sort -u > ./RGs/${SM}.numlibs

done

rm ./RGs/*tmp*
rm ./RGs/*txt*
```

# Identify Neutral Sites

Take the gff file, remove the whole 'region' (chromosome length), the remaining space is intergenic:  

```bash
bedtools sort -i GCA_017976375.1_bCucCan1.pri_genomic.CHR.gff | awk '$3 != "region"' | bedtools complement -i - -g GCA_017976375.1_bCucCan1.pri_genomic.CHR.genome.ordered -L > intergenic.bed
```

Grab introns:

```bash
#extract exons
awk '$3=="exon"' GCA_017976375.1_bCucCan1.pri_genomic.CHR.gff | bedtools sort > exons.bed

#extract genes
awk '$3=="gene"' GCA_017976375.1_bCucCan1.pri_genomic.CHR.gff | bedtools sort > genes.bed

# Get introns by subtracting exons from genes
bedtools subtract -a genes.bed -b exons.bed > introns.bed
awk '{OFS="\t"}{print $1, $4-1, $5}' introns.bed > introns_zerobased.bed
```

Now find degenerate sites using [this](https://github.com/zhangrengang/degeneracy) github to identify 4-fold sites:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

GFF=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.gff
GENOME=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

python /dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/neutral/degeneracy/get_degeneracy.py $GFF $GENOME cuckoo

#convert to bed 
awk '{OFS="\t"}{print $1, $2-1, $2}' cuckoo.codon3-fold4.pos > 4fold.bed
```

And then just merge the beds, and ensure you merge overlapping regions:

```bash
cat intergenic.bed introns_zerobased.bed 4fold.bed | \
	bedtools sort -i - | \
	bedtools merge -i - > /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed
	
awk '{print $3- $2}' neutral-sites_cuckoo__intergenic-intron-4fold.bed | datamash sum 1
992158715
```

## Repeat Masking

Create de novo repeat libraries 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=200:00:00

RUN=$1

#conda activate repmod
#RepeatModeler version 2.0.3
BuildDatabase -name ${RUN} ${RUN}.chr.fa
RepeatModeler -database ${RUN} -pa 20 -LTRStruct >& ${RUN}.out
```

Filter TEs, largely following [this](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-021-00259-7#MOESM1) tutorial:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

#conda activate te_ann
#USES cd-hit and pfam to identify redundant repeats and any which overlap proteins

pfamdb=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/wchrom/repeats/2023mar/pfam

#ensure there's no duplicate sequence IDs
seqkit rename Avian-families.fa > Avian-families.id.fa

#identify redundant hits
cd-hit-est -T 10 -i Avian-families.id.fa -o Avian-cdhit.fa -M 0 -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500

#grab longest ORF
python -m jcvi.formats.fasta longestorf Avian-cdhit.fa > Avian-cdhit.orf.fasta

#translate as if they are vertebrate CDS
transeq -sequence Avian-cdhit.orf.fasta -outseq Avian-cdhit.orf.trns.fa -table 0 -frame 1

#and scan for proteins
pfam_scan.pl -fasta Avian-cdhit.orf.trns.fa -dir $pfamdb > pfam.results
#no hits found 
```

Repeatmask by compartment (Z/W/A):

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=200:00:00

#conda activate repmod
RUN=$1

#Now run repeatmasker with that library against the autosomes, W, and Z independently
for CHR in $(cat CHRS.list); do
seqtk subseq ${RUN}.chr.fa ${CHR}.list > ${RUN}-${CHR}.fa
RepeatMasker -pa 20 -a -s -gff -no_is -lib Avian-cdhit.fa ${RUN}-${CHR}.fa &> RMaves_${RUN}-${CHR}.run.out
done
```

Parse outputs:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=3
#SBATCH --time=24:00:00

#conda activate orthofinder
RUN=$1

#Now run repeatmasker with that library against the autosomes, W, and Z independently
for CHR in $(cat CHRS.list); do
perl /dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/repeatmodeler/Parsing-RepeatMasker-Outputs/parseRM.pl -i ${RUN}-${CHR}.fa.align -p -f ${RUN}.chr.fa -r Avian-cdhit.fa -v
done
```

And now, grab the repeats in bed format:

```bash
awk '{OFS="\t"}{print "chr_"$5, $6, $7, $9, $10, $11}' *txt | bedtools sort -i - > ~/symlinks/c00c00/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats.bed

#but, let's retain low complexity and simple repeats, and allow MQ to deal with that:
egrep -v 'Simple|Low_complexity' GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats.bed > GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed
```





# Spatial

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/miganc')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(sf)
library(ggspatial)
library(spThin)
library(factoextra)
library(ggforce)
library(ggsci)

#read metadata
md = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/miganc/Full_Metadata.txt')

#set up map and convert df to coordinate frame
world = map_data("world")
md = md %>% mutate(LatJit = jitter(Latitude,amount =1),
                    LonJit = jitter(Longitude,amount=1))
sites = st_as_sf(md, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant") 
spp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=Species,shape=Species),
          size=2.5,alpha=0.8,show.legend = T,stroke=0.5) +
  scale_fill_aaas()+
  scale_shape_manual(values=c(21,24))+
  xlab('') + ylab('') +
  coord_sf(xlim = c(min(md$Longitude) - 15, max(md$Longitude) + 5), 
           ylim = c(min(md$Latitude) - 65, max(md$Latitude) + 5), expand = FALSE) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.position = 'top') +
  annotation_scale(line_width = 0.5)

pdf('figures/AllSamples_Extent_2023NOV22.pdf',height=6,width=9)
spp
dev.off()
```

Plot migration data, popgen data:

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/miganc')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(sf)
library(ggspatial)
library(spThin)
library(factoextra)
library(ggforce)
library(ggsci)
library(lubridate)
library(geosphere)
library(concaveman) 

#read metadata
md = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/miganc/Full_Metadata.txt')
md = md %>% mutate(Species = gsub('C. canorus','Cuculus canorus',Species),
                   Species = gsub('C. optatus','Cuculus optatus',Species))

#read spatial data 
m = read_csv('03-BoD-deaths-removed-no-NDVI.csv')
m %>% select(dataset, individual.taxon.canonical.name,tag.local.identifier) %>% unique %>% count(individual.taxon.canonical.name)

#set up map and convert df to coordinate frame
world = map_data("world")
md = md %>% mutate(LatJit = jitter(Latitude,amount =1),
                   LonJit = jitter(Longitude,amount=1))
sites = st_as_sf(md, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant") 

tracks_df = m %>% filter(!grepl('satura',individual.taxon.canonical.name)) %>% select(individual.taxon.canonical.name,timestamp,Route,Wolf,tag.local.identifier,location.long,location.lat)
names(tracks_df) = c('Species','DT','Route','WolfID','ID','Long','Lat')

#assign 'year1', 'year2', 'year3', etc 
tracks_df = tracks_df %>%
  group_by(ID) %>%
  mutate(DT = mdy_hm(DT)) %>% 
  arrange(DT) %>%
  mutate(
    Year_Num = cumsum(month(DT) == 11 & lag(month(DT), default = 10) != 11),  # Increment each November
    Year_Label = paste0("Year", Year_Num)
  ) %>%
  ungroup()

tracks_df = tracks_df %>%
  mutate(
    Status = case_when(
      month(DT) %in% c(12) ~ "Overwintering",
      month(DT) %in% c(6) ~ "Breeding",
      TRUE ~ NA_character_
    )
  )


# #remove any points further than 1000km 
# tracks_df <- tracks_df %>%
#   group_by(ID) %>%
#   arrange(DT) %>%
#   mutate(
#     Prev_Long = lag(Long),
#     Prev_Lat = lag(Lat),
#     Distance = distHaversine(cbind(Prev_Long, Prev_Lat), cbind(Long, Lat))
#   ) %>%
#   ungroup() %>% 
#   mutate(Distance = Distance / 1000)
# 
# tfx = tracks_df %>% drop_na(Distance) %>% 
#   filter(Distance < 250) #1000km x 1000m

#create points 
tracks = st_as_sf(tracks_df, coords = c("Long", "Lat"), 
                 crs = 4326, agr = "constant") 
#create lines 
tracks_lines <- tracks %>% 
  group_by(ID, Species) %>% 
  filter(n() >= 2) %>%  # Keep only groups with at least two records
  summarize(geometry = st_union(geometry)) %>%
  st_cast("LINESTRING") %>%
  st_as_sf()

#create convex hull polygons 
#could only install concaveman and V8 with this:
#Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
#install.packages('V8')
counter = 0 
for (stat in c('Breeding','Overwintering')) { for (sp in c('Cuculus canorus','Cuculus optatus')) { counter = counter + 1
  t = tracks %>% filter(Status == stat & Species == sp) 
  p = concaveman(t,concavity=0.5,length_threshold=25)
  p = p %>% mutate(Species = sp, Status = stat)
  assign(paste0('p',counter),p)
}}
ps = rbind(p1,p2,p3,p4)

tpp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data=tracks_lines,col='grey20',alpha=0.2,
          size=0.15)+
  # geom_sf(data = tracks %>% filter(Status == 'Breeding' | Status == 'Overwintering'),
  #         aes(fill=Status,shape=Species),
  #         size=2.5,alpha=0.8,show.legend = T,stroke=0.5) +
  geom_sf(data = ps,aes(fill=Status),alpha=0.5)+
  geom_sf(data = sites,pch=21,fill='white',col='black',size=2.5,alpha=0.8,show.legend = T,stroke=0.5) +
  scale_fill_aaas()+
  scale_color_brewer(palette = 'Set2')+
  scale_shape_manual(values=c(16,15))+
  facet_wrap(Species~.)+
  xlab('') + ylab('') +
  coord_sf(xlim = c(min(md$Longitude) - 15, max(md$Longitude) + 5), 
           ylim = c(min(md$Latitude) - 65, max(md$Latitude) + 5), expand = FALSE) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = 'top') +
  annotation_scale(line_width = 0.5)
tpp
pdf('figures/AllSamples_Extent-Species_2023NOV27.pdf',height=6,width=12)
tpp
dev.off()

#only BTO
uk_filtered <- tracks_df %>%
  filter(Status == "Breeding", 
         Lat >= 49, Lat <= 61,   # UK latitude boundaries
         Long >= -8, Long <= 2)  # UK longitude boundaries

#filter on distance 
tracks_df <- tracks_df %>%
  group_by(ID) %>%
  arrange(DT) %>%
  mutate(
    Prev_Long = lag(Long),
    Prev_Lat = lag(Lat),
    Distance = distHaversine(cbind(Prev_Long, Prev_Lat), cbind(Long, Lat)) / 1000,  # Distance in kilometers
    Long_Segment = Distance > 2000  # Identify long segments
  ) %>%
  ungroup() %>% drop_na(Distance)

#flag 
tracks_df <- tracks_df %>%
  group_by(ID) %>%
  mutate(
    Remove_Point = lag(Long_Segment, default = FALSE) | Long_Segment | lead(Long_Segment, default = FALSE)
  ) %>%
  ungroup()

tfx = tracks_df %>%
  filter(!Remove_Point)

bto_tracks = tracks_df %>% filter(ID %in% uk_filtered$ID)
#label as spring or fall migration
bto_tracks = bto_tracks %>%
  mutate(Status = case_when(
    month(DT) %in% c(3,4,5) ~ "Spring",
    month(DT) %in% c(7,8,9) ~ "Fall",
    TRUE ~ NA))
btosf =  st_as_sf(bto_tracks, coords = c("Long", "Lat"), 
                  crs = 4326, agr = "constant") 
bto_lines <- btosf %>% 
  group_by(ID,Status,Route,WolfID) %>% 
  filter(n() >= 2) %>%  # Keep only groups with at least two records
  summarize(geometry = st_union(geometry)) %>%
  st_cast("LINESTRING") %>%
  st_as_sf()
bto_mig = bto_lines %>% filter(!is.na(Status)) %>% filter(Route == 'West' | Route == 'East') %>% filter(!is.na(WolfID))
bbox <- st_bbox(btosf)
up = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  #geom_sf(data=bto_tracks,col='grey20',alpha=0.2,size=0.5)+
  geom_sf(data=bto_mig,aes(col=Route),alpha=0.5,size=0.15)+
  xlab('') + ylab('') +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE)+
  theme_classic() +
  facet_wrap(Status~.)+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = 'top') +
  annotation_scale(line_width = 0.5)

pdf('figures/UKSamples_Route_2023NOV27.pdf',height=6,width=12)
up
dev.off()

```

## Plot Spherical

Plotting on a spherical globe is easy using this [code](https://stackoverflow.com/questions/70756215/plot-geodata-on-the-globe-perspective-in-r). 

```R
#### Plot Track data on 3D Globe 
setwd('~/merondun/cuculus_migration/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(sf)
library(tidyverse)
library(ggspatial)
library(RColorBrewer)
library(giscoR)
library(ggpubr)

# Read metadata
md <- read_tsv('Full_Metadata.txt') %>%
  mutate(
    Species = str_replace_all(Species, c('C. canorus' = 'Cuculus canorus', 'C. optatus' = 'Cuculus optatus')),
    LatJit = jitter(Latitude, amount = 1),
    LonJit = jitter(Longitude, amount = 1)
  )

# Read spatial data, filter according to Maire's filtering
final.locations <- read_csv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/spatial/03-BoD-deaths-removed-no-NDVI.csv')

# Calculate distance between each point
final.locations <- final.locations %>%
  group_by(tag.local.identifier) %>%
  mutate(long.to.sf = location.long,
         lat.to.sf = location.lat) %>%
  st_as_sf(coords = c("long.to.sf", "lat.to.sf"), crs = 4326) %>%
  arrange(tag.local.identifier, as.POSIXct(timestamp, format = "%Y-%m-%d")) %>%
  mutate(dist = c(0,as.numeric(st_distance(geometry[-1], geometry[-n()], by_element = T))/1000)) %>%
  st_drop_geometry()

# Remove paths over 3000km?
large.gap <- which(final.locations$dist > 3000) # Which paths are large?
final.locations$large.gap <- F # Identify previous point before large path as cut off
final.locations$large.gap[large.gap-1] <- T

# Create function to identify large gap
gap.function <- function(var){
  x <- rle(var)
  cumsum(x[[2]]) |> (\(.){ .[which(!x[[2]], T)] <- NA ; .})() |> rep(x[[1]])
}

# Create new ID when gap is over 3000
final.locations <- final.locations %>%
  group_by(tag.local.identifier) %>%
  mutate(new_grp = gap.function(large.gap)) %>%
  fill(new_grp, .direction = "up") %>%
  mutate(new_id = paste0(tag.local.identifier, new_grp))

# Plot paths
ggplot() +
  geom_path(final.locations, mapping = aes(location.long, location.lat, group = new_id))

# Create back-up, switch to Justin workflow and rename sensible columns
m <- final.locations
names(m)[names(m) == 'individual.taxon.canonical.name'] <- 'Species'
names(m)[names(m) == 'new_id'] <- 'ID'

# Add R interpretable datestring, add 'Breeding/Overwintering
tracks_df <- m %>%
  filter(!str_detect(Species, 'satura')) %>%
  select(Species, timestamp, Route, Wolf, ID, location.long, location.lat) %>%
  rename(
    DT = timestamp,
    Long = location.long,
    Lat = location.lat
  ) %>%
  mutate(
    DT = mdy_hm(DT),
    Year_Num = cumsum(month(DT) == 11 & lag(month(DT), default = 10) != 11),
    Year_Label = paste0("Year", Year_Num),
    Status = case_when(
      month(DT) %in% c(12) ~ "Overwintering",
      month(DT) %in% c(6) ~ "Breeding",
      TRUE ~ NA_character_
    )
  )

# Assign CCW / CCE / COW / COE based on longitude during June and July
# routes <- tracks_df %>%
#   filter(month(DT) == 6 | month(DT) == 7) %>% 
#   group_by(tag.local.identifier,Species) %>% 
#   summarize(Long = mean(Long), Lat = mean(Lat)) %>% 
#   mutate(Direction = ifelse( Species == 'Cuculus canorus' & Long > 110, 'CCE',
#                              ifelse( Species == 'Cuculus canorus','CCW',
#                                      ifelse( Species == 'Cuculus optatus' & Long > 115,'COE','COW')))) %>% 
#   select(-c(Long,Lat))

# Assign CCW / CCE / CO - 3 pop 
routes <- tracks_df %>%
  filter(month(DT) == 6 | month(DT) == 7) %>% 
  group_by(tag.local.identifier,Species) %>% 
  summarize(Long = mean(Long), Lat = mean(Lat)) %>% 
  mutate(Direction = ifelse( Species == 'Cuculus canorus' & Long > 110, 'CCE',
                             ifelse( Species == 'Cuculus canorus','CCW', 'CO')))%>% 
  select(-c(Long,Lat))

# BUT, some indivduals dont' have track data in June July, manually inspect and assign - they are all western
tracks_routes <- left_join(tracks_df, routes) %>% 
  mutate(Direction = ifelse(!is.na(Direction),Direction,
                            ifelse(Species == 'Cuculus canorus','CCW','CO')))

# Confirm 
tracks_routes %>% 
  ggplot() +
  geom_path(mapping = aes(Long, Lat, group = ID, col = Direction)) +
  facet_grid(Species~.,scales='free')

# Convert to sf
tracks = st_as_sf(tracks_routes, coords = c("Long", "Lat"), 
                  crs = 4326, agr = "constant") 

# Colors for can and opt 
cols = brewer.pal(8,'Paired')[c(3,4,7,8)]

# CANORUS: projection string used for the polygons & ocean background
crs_string <- "+proj=ortho +lon_0=60 +lat_0=30"

# background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# country polygons, cut to size
world <- gisco_countries %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) # reproject to ortho


# Reproject into view of africa / asia 
reproj_can <- tracks %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) %>% # reproject to ortho 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2]) %>% 
  filter(Species == 'Cuculus canorus') %>% 
  st_transform(crs = crs_string)

# now the action!
can_plot <- ggplot(data = world) +
  geom_sf(data = ocean, fill = "cadetblue1", color = NA) + # background first
  geom_sf(lwd = .1,col='grey80',fill='grey99') + # now land over the oceans
  geom_path(reproj_can, mapping = aes(x = Long, y = Lat, group = ID, col = Direction),
            alpha=0.5) + #and the cuckoos
  scale_color_manual(values=cols[1:2])+
  theme_void()
can_plot

# OPTATUS: projection string used for the polygons & ocean background
crs_string <- "+proj=ortho +lon_0=100 +lat_0=30"

# background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# country polygons, cut to size
world <- gisco_countries %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) # reproject to ortho

# Reproject into view of africa / asia 
reproj_opt <- tracks %>% 
  filter(Species == 'Cuculus optatus') %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) %>% # reproject to ortho 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2]) %>% 
  st_transform(crs = crs_string)

# now the action! 
opt_plot <- ggplot(data = world) +
  geom_sf(data = ocean, fill = "cadetblue1", color = NA) + # background first
  geom_sf(lwd = .1,col='grey80',fill='grey99') + # now land over the oceans
  geom_path(reproj_opt, mapping = aes(x = Long, y = Lat, group = ID, col = Direction),
            alpha=0.8) + # More opaque for optatus to help distinguish
  scale_color_manual(values=cols[4])+  #and the colors
  theme_void()
opt_plot

ggsave('figures/20241218_Spherical-Tracks-Filtered-CanOpt-3Pop.pdf',
       ggarrange(can_plot,opt_plot),height=4,width=9,dpi=600)

saveRDS(tracks, "/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/spatial/tracks_data.rds")
```

![image-20241022163640949](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20241022163640949.png)

### Spherical Satellite

```bash
#### Plot Track data on 3D Globe 
setwd('~/merondun/cuculus_migration/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(qdrmapbox)
library(raster)
library(sf)
library(tidyverse)
library(ggspatial)
library(RColorBrewer)
library(giscoR)
library(ggpubr)

# Read metadata
md <- read_tsv('Full_Metadata.txt') %>%
  mutate(
    Species = str_replace_all(Species, c('C. canorus' = 'Cuculus canorus', 'C. optatus' = 'Cuculus optatus')),
    LatJit = jitter(Latitude, amount = 1),
    LonJit = jitter(Longitude, amount = 1)
  )

# Read spatial data, filter according to Maire's filtering
final.locations <- read_csv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/spatial/03-BoD-deaths-removed-no-NDVI.csv')

# Calculate distance between each point
final.locations <- final.locations %>%
  group_by(tag.local.identifier) %>%
  mutate(long.to.sf = location.long,
         lat.to.sf = location.lat) %>%
  st_as_sf(coords = c("long.to.sf", "lat.to.sf"), crs = 4326) %>%
  arrange(tag.local.identifier, as.POSIXct(timestamp, format = "%Y-%m-%d")) %>%
  mutate(dist = c(0,as.numeric(st_distance(geometry[-1], geometry[-n()], by_element = T))/1000)) %>%
  st_drop_geometry()

# Remove paths over 3000km?
large.gap <- which(final.locations$dist > 3000) # Which paths are large?
final.locations$large.gap <- F # Identify previous point before large path as cut off
final.locations$large.gap[large.gap-1] <- T

# Create function to identify large gap
gap.function <- function(var){
  x <- rle(var)
  cumsum(x[[2]]) |> (\(.){ .[which(!x[[2]], T)] <- NA ; .})() |> rep(x[[1]])
}

# Create new ID when gap is over 3000
final.locations <- final.locations %>%
  group_by(tag.local.identifier) %>%
  mutate(new_grp = gap.function(large.gap)) %>%
  fill(new_grp, .direction = "up") %>%
  mutate(new_id = paste0(tag.local.identifier, new_grp))

# Plot paths
ggplot() +
  geom_path(final.locations, mapping = aes(location.long, location.lat, group = new_id))

# Create back-up, switch to Justin workflow and rename sensible columns
m <- final.locations
names(m)[names(m) == 'individual.taxon.canonical.name'] <- 'Species'
names(m)[names(m) == 'new_id'] <- 'ID'

# Add R interpretable datestring, add 'Breeding/Overwintering
tracks_df <- m %>%
  filter(!str_detect(Species, 'satura')) %>%
  select(Species, timestamp, Route, Wolf, ID, location.long, location.lat) %>%
  rename(
    DT = timestamp,
    Long = location.long,
    Lat = location.lat
  ) %>%
  mutate(
    DT = mdy_hm(DT),
    Year_Num = cumsum(month(DT) == 11 & lag(month(DT), default = 10) != 11),
    Year_Label = paste0("Year", Year_Num),
    Status = case_when(
      month(DT) %in% c(12) ~ "Overwintering",
      month(DT) %in% c(6) ~ "Breeding",
      TRUE ~ NA_character_
    )
  )

# Assign CCW / CCE / COW / COE based on longitude during June and July
# routes <- tracks_df %>%
#   filter(month(DT) == 6 | month(DT) == 7) %>% 
#   group_by(tag.local.identifier,Species) %>% 
#   summarize(Long = mean(Long), Lat = mean(Lat)) %>% 
#   mutate(Direction = ifelse( Species == 'Cuculus canorus' & Long > 110, 'CCE',
#                              ifelse( Species == 'Cuculus canorus','CCW',
#                                      ifelse( Species == 'Cuculus optatus' & Long > 115,'COE','COW')))) %>% 
#   select(-c(Long,Lat))

# Assign CCW / CCE / CO - 3 pop 
routes <- tracks_df %>%
  filter(month(DT) == 6 | month(DT) == 7) %>% 
  group_by(tag.local.identifier,Species) %>% 
  summarize(Long = mean(Long), Lat = mean(Lat)) %>% 
  mutate(Direction = ifelse( Species == 'Cuculus canorus' & Long > 110, 'CCE',
                             ifelse( Species == 'Cuculus canorus','CCW', 'CO')))%>% 
  select(-c(Long,Lat))

# BUT, some indivduals dont' have track data in June July, manually inspect and assign - they are all western
tracks_routes <- left_join(tracks_df, routes) %>% 
  mutate(Direction = ifelse(!is.na(Direction),Direction,
                            ifelse(Species == 'Cuculus canorus','CCW','CO')))

# Confirm 
tracks_routes %>% 
  ggplot() +
  geom_path(mapping = aes(Long, Lat, group = ID, col = Direction)) +
  facet_grid(Species~.,scales='free')

# Convert to sf
tracks = st_as_sf(tracks_routes, coords = c("Long", "Lat"), 
                  crs = 4326, agr = "constant") 

# Colors for can and opt 
cols = brewer.pal(8,'Paired')[c(3,4,7,8)]

# CANORUS: projection string used for the polygons & ocean background
crs_string <- "+proj=ortho +lon_0=60 +lat_0=30"

# background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# country polygons, cut to size
world <- gisco_countries %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) # reproject to ortho


# Reproject into view of africa / asia 
reproj_can <- tracks %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) %>% # reproject to ortho 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2]) %>% 
  filter(Species == 'Cuculus canorus') %>% 
  st_transform(crs = crs_string)

# now the action!
can_plot <- ggplot(data = world) +
  geom_sf(data = ocean, fill = "cadetblue1", color = NA) + # background first
  geom_sf(lwd = .1,col='grey80',fill='grey99') + # now land over the oceans
  geom_path(reproj_can, mapping = aes(x = Long, y = Lat, group = ID, col = Direction),
            alpha=0.5) + #and the cuckoos
  scale_color_manual(values=cols[1:2])+
  theme_void()
can_plot

# OPTATUS: projection string used for the polygons & ocean background
crs_string <- "+proj=ortho +lon_0=100 +lat_0=30"

# background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# country polygons, cut to size
world <- gisco_countries %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) # reproject to ortho

# Reproject into view of africa / asia 
reproj_opt <- tracks %>% 
  filter(Species == 'Cuculus optatus') %>% 
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = crs_string) %>% # reproject to ortho 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2]) %>% 
  st_transform(crs = crs_string)

# now the action! 
opt_plot <- ggplot(data = world) +
  geom_sf(data = ocean, fill = "cadetblue1", color = NA) + # background first
  geom_sf(lwd = .1,col='grey80',fill='grey99') + # now land over the oceans
  geom_path(reproj_opt, mapping = aes(x = Long, y = Lat, group = ID, col = Direction),
            alpha=0.8) + # More opaque for optatus to help distinguish
  scale_color_manual(values=cols[4])+  #and the colors
  theme_void()
opt_plot

ggsave('figures/20241218_Spherical-Tracks-Filtered-CanOpt-3Pop.pdf',
       ggarrange(can_plot,opt_plot),height=4,width=9,dpi=600)

saveRDS(tracks, "/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/spatial/tracks_data.rds")


### Satelltie
library(ggmap)
register_google(key = "AIzaSyDzjCgFf9hhFIBkVJ4tAyy6cjjq9BjmG4U")

# Canorus 
can_flat <- tracks %>% 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2]) %>% 
  filter(Species == 'Cuculus canorus') 

can_bb <- st_bbox(can_flat)
names(can_bb) <- c('left','bottom','right','top')
can_bb_buffered <- can_bb + c(-buffer_sides, -buffer_tops, buffer_sides, buffer_tops)

map_data <- get_map(location = can_bb, zoom = 2, source = "google", maptype = 'satellite')
base_map <- ggmap(map_data)
base_map

# now the action!
can_plot <- base_map +
  borders("world", col = "gray90", size = 0.05,alpha=0.5)+
  geom_path(can_flat, mapping = aes(x = Long, y = Lat, group = ID, col = Direction),
            alpha=0.5) + #and the cuckoos
  scale_color_manual(values=cols[1:2])+
  theme_void()
can_plot


# Optatus
opt_flat <- tracks %>% 
  mutate(Long = st_coordinates(geometry)[,1],
         Lat = st_coordinates(geometry)[,2]) %>% 
  filter(Species == 'Cuculus optatus') 

opt_bb <- st_bbox(opt_flat)
names(opt_bb) <- c('left','bottom','right','top')

map_data <- get_map(location = opt_bb, zoom = 2, source = "google", maptype = 'satellite')
base_map <- ggmap(map_data)
base_map

# now the action!
opt_plot <- base_map +
  borders("world", col = "gray90", size = 0.05,alpha=0.5)+
  geom_path(opt_flat, mapping = aes(x = Long, y = Lat, group = ID, col = Direction),
            alpha=0.5,col=cols[4]) + #and the cuckoos
  theme_void()
opt_plot

both <- ggarrange(can_plot,opt_plot)
ggsave('figures/20241218_Satellite_3Pop.pdf',dpi=300,both,height=6,width=9)

```

# Plot Bioacoustic

Plot syllable duration and pitch:

```bash
setwd('C:/Users/herit/Dropbox/CUCKOO_migration/Manuscript/Figures/Figure_1_TracksPopGenStructure/')
library(tidyverse)
library(RColorBrewer)

sounds <- read_csv('Cuckoos sounds/sounds.csv')

# Retain species
analyze <- sounds %>% filter(grepl('can|opt',Species))
cols <- c('#b2df8a','#33a02c','#fe9024')

pca_plot <- analyze %>% 
  ggplot() +
  geom_point(aes(x=`Dur of syllable`,y=`Pitch of syll`,fill = Species,shape=Species), 
             size = 2, color = 'grey20', alpha = 0.9) +
  theme_bw(base_size=8) +
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,21,24))+
  xlab(paste0('Duration of Syllable (s)'))+
  ylab(paste0('Pitch of Syllable (Hz)'))
pca_plot

ggsave('Cuckoos sounds/20250212_Simple2Var-Sounds-CanOpt.pdf',units='in',dpi=300,height=2,width=3)

```

PCA:

```R
setwd('C:/Users/herit/Dropbox/CUCKOO_migration/Manuscript/Figures/Figure_1_TracksPopGenStructure/')
library(tidyverse)
library(RColorBrewer)

scores <- read_csv('Cuckoos sounds/PCA_Scores.csv')
ids <- read_csv('Cuckoos sounds/sounds.csv')
vals <- read_csv('Cuckoos sounds/importance 2s.csv')
vt <- vals %>% filter(grepl('Proportion of',`...1`)) %>% select(-`...1`) %>% t() %>% as.data.frame
vt$Axis <- rownames(vt)
rownames(vt) <- NULL

p <- left_join(scores %>% dplyr::rename(Number = `...1`),ids)
cols <- c('#b2df8a','#33a02c','#fe9024')

pca_plot <- p %>% 
  ggplot() +
  geom_point(aes(x=PC1,y=-PC2,fill = Species,shape=Species), 
             size = 2, color = 'grey20', alpha = 0.9) +
  theme_bw(base_size=8) +
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,21,24))+
  xlab(paste0('PC1 (',round(vt[1,1],3)*100,'%)'))+
  ylab(paste0('PC2 (',round(vt[2,1],3)*100,'%)'))
pca_plot

ggsave('Cuckoos sounds/20250115_PCA-Sounds-CanOpt.pdf',units='in',dpi=300,height=2,width=3)
```

# Plot Hindcasting

```R
setwd('C:/Users/herit/Dropbox/CUCKOO_migration/Manuscript/Figures/Figure_3_Hindcasting/Files for Fig 3/')
library(tidyverse)
library(readxl)
library(ggforce)
library(RColorBrewer)

# Load the Excel file, skipping the first row
df <- read_excel("BreedWinterOverlap3.xlsx",sheet = 'R_Input') %>% as_tibble()
df <- df %>% mutate(year_ka = year / 1000)

# Add temp periods (For demography! )
periods <- tibble(
  xmin = c(0, 23750, 121500, 56750),
  xend = c(10000, 33750, 131500, 66750),
  col = c("hot", "cold", "hot", "cold"),
  label = c("t0", "t0", "t1", "t1"),
  y = c(0.005, -7.6, 1.94, -6.46)
) %>% 
  mutate(xmin = xmin / 1000,
         xend = xend / 1000)

# Temperature 
temp <- ggplot(df, aes(x = year_ka, y = Temp, color = Temp)) +
  geom_line() +
  scale_color_gradientn(colors = c("blue", "red")) +
  geom_rect(data=periods,aes(xmin=xmin,xmax=xend,ymin=-Inf,ymax=Inf,fill=col),inherit.aes=FALSE,alpha=0.25)+
  scale_fill_manual(values=c('cyan3','salmon2'))+
  scale_x_reverse(limits=c(800,0),
                  breaks = seq(0, 800, by = 100),  
                  labels = function(x) paste0(x, " Ka")) +
  labs(y = "Temperature", x = "Year") +
  facet_grid(1~.,scales='free')+
  theme_bw(base_size=8) +
  theme(legend.position='none')
temp

ggsave('20250206_Temperature.pdf',temp,height=0.65,width=6.5,dpi=600,units='in')

# Gap between breeding / overwintering
df_long <- df %>% dplyr::rename(East = east_diff, West = west_diff) %>% 
  pivot_longer(
    cols = c(West, East),
    names_to = "region",
    values_to = "difference"
  )

breeding <- ggplot(df_long, aes(x = year_ka, y = difference, group = region)) +  
  geom_link2(aes(colour = after_stat(ifelse(y > 0, "Overlap", "Gap"))),
             linewidth = 0.5) +
  scale_x_reverse(limits=c(800,0),
    breaks = seq(0, 800, by = 100),  
    labels = function(x) paste0(x, " Ka")) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_rect(data=periods,aes(xmin=xmin,xmax=xend,ymin=-Inf,ymax=Inf,fill=col),inherit.aes=FALSE,alpha=0.25)+
  scale_fill_manual(values=c('cyan3','salmon2'))+
  scale_color_manual(values=c('darkorchid','seagreen','darkblue','yellow'))+
  labs(y = "Overlap (+) or Gap (-)", x = "Year") +
  scale_y_continuous(limits = c(-40, 40)) +
  #facet_grid(region ~ ., scales = "free") + #for separated east/west
  theme_bw(base_size = 8)+
  theme(legend.position='none')
breeding

df_long <- df_long %>%
  mutate(color = case_when(
    region == "East" & difference >= 0 ~ "A",
    region == "West" & difference >= 0 ~ "B",
    region == "East" & difference < 0 ~ "C",
    region == "West" & difference < 0 ~ "D"
  ))

breeding <- ggplot(df_long, aes(x = year_ka, y = difference, group = region)) +  
  geom_link2(aes(colour = color), linewidth = 0.5) +
  scale_x_reverse(limits = c(800, 0),
                  breaks = seq(0, 800, by = 100),  
                  labels = function(x) paste0(x, " Ka")) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_rect(data = periods, aes(xmin = xmin, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), inherit.aes = FALSE, alpha = 0.25) +
  scale_fill_manual(values = c('cyan3', 'salmon2')) +
  scale_color_manual(values = c('darkblue', 'blue', 'darkgreen', 'lightgreen')) +
  labs(y = "Overlap (+) or Gap (-)", x = "Year") +
  scale_y_continuous(limits = c(-40, 40)) +
  theme_bw(base_size = 8) +
  theme(legend.position = 'none')
breeding

# The line itself is difficult to vectorize because it is thousands of segments, export once as png, and then hash out the link() lines and resave
ggsave('20250206_BreedingOverlap.png',breeding,height=0.95,width=6.5,dpi=600,units='in')
ggsave('20250206_BreedingOverlapAxes.pdf',breeding,height=0.95,width=6.5,dpi=600,units='in')


# Gaps
df_tracks <- df %>% 
  pivot_longer(c(EWgap, Sahara0Gap,IndianOcean)) 
df_tracks$name <- factor(df_tracks$name,levels=c('EWgap','Sahara0Gap','IndianOcean'))  
tracks <- df_tracks %>% 
  ggplot(aes(x = year_ka, y = value,col=name)) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept=0.5,lty=2,col='grey30')+
  scale_x_reverse(limits=c(800,0),
                  breaks = seq(0, 800, by = 100),  
                  labels = function(x) paste0(x, " Ka")) +
  geom_rect(data=periods,aes(xmin=xmin,xmax=xend,ymin=-Inf,ymax=Inf,fill=col),inherit.aes=FALSE,alpha=0.25)+
  scale_fill_manual(values=c('cyan3','salmon2'))+
  labs(y = "Temperature", x = "Year")+
  scale_color_manual(values=brewer.pal(3,'Set2'))+
  #facet_grid(name~.,scales='free')+
  theme_bw(base_size=8) +
  coord_cartesian(ylim=c(0,1))
  theme(legend.position='none')
tracks

ggsave('20250206_Gaps.pdf',tracks,height=0.95,width=6.5,dpi=600,units='in')


```

# Temperature Periods 

Identify warm and cold periods: 

   ```
.libPaths('~/mambaforge/envs/R/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_chyiyin/crosscoal/output')
library(tidyverse)

#import ice core data 
icecore = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Ice_Temperature_Reconstructions_Kawamura_NOAA-6076.txt')
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
                                                label=paste0(stage,': ',period,'  ',round(temp,2),'°C\n',time_low/1000,' - ',time_high/1000,'Ka')),col='black')+
  xlab('Time')+ylab('Delta T (°C)')+
  theme_test(base_size=10)
ice_plot

png('~/merondun/cuculus_migration/figures/WarmCold_Periods.png',units='in',res=300,height=2,width=3.5)
pdf('~/merondun/cuculus_migration/figures/WarmCold_Periods.pdf',height=2,width=3.5)
ice_plot
dev.off()

   ```



| Start  | End    | Duration | Stage | Period | Mean Temperature |
| ------ | ------ | -------- | ----- | ------ | ---------------- |
| 0      | 10000  | 10000    | hot   | t0     | 0.005            |
| 23750  | 33750  | 10000    | cold  | t0     | -7.6             |
| 121500 | 131500 | 10000    | hot   | t1     | 1.94             |
| 56750  | 66750  | 10000    | cold  | t1     | -6.46            |

​     

   

# Processing

## Alignment

Still on the 658 individual SRR accessions. Submit positional library.

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

SCRATCH=tmp/$SLURM_JOB_ID

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq
qcdata=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/2021_04
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

RUN=$1

bwa index ${genome}
#Mapping script, will map reads with the prefix identified in the list ${lists}/${RUN}

        ID="$(cat ${wd}/RGs/${RUN}.ID | tail -1)";
        PU="$(cat ${wd}/RGs/${RUN}.PU | tail -1)";
        SM="$(cat ${wd}/RGs/${RUN}.SM | tail -1)";
        LB="$(cat ${wd}/RGs/${RUN}.LB | tail -1)";

        #library mapping
        bwa mem -M -p -t 10 -R "@RG\tID:${ID}\tSM:${SM}\tPL:ILLUMINA\tPU:${PU}\tLB:${LB}" ${genome} ${qcdata}/${RUN}.trim.fastq.gz | samtools sort -@10 -o ${SCRATCH}/${RUN}.bam -;
        samtools view -b -f 1 -F 524 ${SCRATCH}/${RUN}.bam > ${outdir}/${RUN}.bam
        samtools index -b ${outdir}/${RUN}.bam;
```

And submit by sample:

```bash
for sample in $(cat Libraries.list); do sbatch -J ${sample} Align.sh ${sample}; done 
```

Plot Alignment:

```R
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

```


## Merge by Sample

This script will merge, mark duplicates, remove reads overhanging scaffold ends, and summarize (coverage + alignments). Submit positional SAMPLE (n=390). 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

SCRATCH=tmp/$SLURM_JOB_ID

#Bam cleaning script. Merge by sample, Mark duplicates, remove reads overhanging scaffolds, and summarize

RUN=$1

RGdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/trimmed_fastq/RGs
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/2021_04
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir tmp
mkdir $SCRATCH

#Make final bam and stats folder
mkdir ${bamdir}/merged
mkdir ${bamdir}/merged/stats
mkdir ${bamdir}/merged/stats/alignments
mkdir ${bamdir}/merged/stats/coverage

#merge all the libraries together
samtools merge - ${bamdir}/${RUN}*.bam | samtools sort -@3 -o $SCRATCH/${RUN}.bam
samtools index $SCRATCH/${RUN}.bam

sambamba markdup --tmpdir $SCRATCH/ $SCRATCH/${RUN}.bam $SCRATCH/${RUN}.dup.bam
samtools index $SCRATCH/${RUN}.dup.bam

#Remove reads over scaffold end
gatk CleanSam -I $SCRATCH/${RUN}.dup.bam -O ${bamdir}/merged/${RUN}.bam -R ${genome} --CREATE_INDEX true

Summary Stats
gatk CollectAlignmentSummaryMetrics -R ${genome} -I ${bamdir}/merged/${RUN}.bam -O ${bamdir}/merged/stats/alignments/${RUN}.alignments.txt -LEVEL SAMPLE -LEVEL READ_GROUP

#calculate coverage
mosdepth --threads 3 --use-median --by 100000 --fast-mode --no-per-base ${bamdir}/merged/stats/coverage/${RUN} ${bamdir}/merged/${RUN}.bam

```

**Sample Stats**

**Alignments**

Check out alignment stats within /stats/alignments/

```bash
mkdir output
grep 'TOTAL' 001_CB_ORW_CHN_F.alignments.txt > output/Alignments.txt

for i in $(ls *.alignments.txt | sed 's/\..*//g'); do 
awk 'NR==13' ${i}.alignments.txt >> output/Alignments.txt
done
```

**Coverage**

Within /coverage/

```bash
#take the unzipped bed.gz file from mosdepth, and change it into the below
for i in $(ls *.regions.bed | rev | cut -c13- | rev); do awk -v s=${i} 'BEGIN{print s}; {print $4}' ${i}.regions.bed > ${i}.cov; done

#create header file with chr / start / end, since all files on the same genome, they all have the same length and this is consistent across samples 
awk '{OFS="\t"}BEGIN{print "CHR","START","END"};{print $1, $2, $3}' 053_CC_RED_FIN_F.regions.bed | sed 's/\ /\t/g' > 0001Header.cov

#combine them all together..
paste *.cov > Master.coverage
```

Also print the GW mean of median:

```bash
for i in $(ls *.regions.bed | rev | cut -c13- | rev); do awk -v s=${i} '{ total += $4 } END { print s, total/NR }' ${i}.regions.bed >> gw_median.coverage ; done
```

And by chromosome:

```bash
for i in $(ls *.summary.txt | sed 's/\..*//g'); do
grep -i -w 'CM030676.1\|CM030677.1\|total\|CM030679.1' ${i}.mosdepth.summary.txt | awk -v s=${i} 'BEGIN{OFS="\t"}{print s, $1, $4}' >> summarize/Compartment.coverage.txt ; done
```

# Variant Calling

Split Genome

```bash
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

#create dict
gatk CreateSequenceDictionary -R $genome
#split genome by Ns 
gatk ScatterIntervalsByNs -R $genome -O scatter.interval_list -OT ACGT

#split intervals
gatk SplitIntervals -R $genome -L scatter.interval_list --scatter-count 100 -O interval_files/
```

## Bcftools

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=240:00:00

vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/2022_11/raw_vcf
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
ints=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/2022_11/snpcalling/intervals

CHR=$1
GROUP=$2

mkdir $vcfdir
if [[ $CHR == 0102 ]]
then
	maxdp=20000
else
	maxdp=100
fi
echo "WORKING ON INTERVAL ${CHR} WITH MAX DEPTH ${maxdp}"
#bcftools 1.16
bcftools mpileup --max-depth ${maxdp} -C 50 -threads 5 -f ${genome} -b ${GROUP}.bam -R $ints/$CHR.bed -a "AD,ADF,ADR,DP,SP" | \
	bcftools call --ploidy 2 --threads 5 -m -Oz -o ${vcfdir}/${CHR}_bcft.vcf.gz
```

Merge by chromosome

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

mkdir ../merged
CHR=$1
echo "${CHR}"
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/${CHR}.consmask
bcftools concat --file-list VCFs.list --threads 5 -a -r ${CHR} -Ou | \
        bcftools sort -Oz -o ../merged/${CHR}.vcf.gz -
bcftools index --threads 5 ../merged/${CHR}.vcf.gz
```

## Assign Distance Groups

Assign samples within 500km of one another as a potential distance group, we will only examine relatedness within these groups.

```R
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
```

## SNP Filtering

Filter at 4 levels. Have this python script available:

```python
import sys
import pysam

# Take inputs from command line arguments
outgroup_samples_file = sys.argv[1]
input_vcf_file = sys.argv[2]
output_prefix = sys.argv[3]

# Load outgroup samples from the input file
with open(outgroup_samples_file, 'r') as file:
    outgroup_samples = [line.strip() for line in file]

# Open the input VCF file
with pysam.VariantFile(input_vcf_file) as vcf_in:
    # Add 'AA' to the INFO fields of the header
    vcf_in.header.info.add('AA', '1', 'String', 'Ancestral Allele')
    # Create an output file
    with pysam.VariantFile(f"{output_prefix}.vcf.gz", 'w', header=vcf_in.header) as vcf_out:
        # Iterate over all records (SNPs)
        for record in vcf_in:
            alleles = []
            for sample in outgroup_samples:
                genotype = record.samples[sample]['GT']
                if genotype.count(None) > 1 or set(genotype) == {0, 1}:
                    # If the genotype is missing or heterozygous, mark this record as unassignable and break
                    alleles = []
                    break
                allele = record.alleles[genotype[0]]
                alleles.append(allele)
            if len(set(alleles)) > 1:
                # If more than one allele is found among the outgroup samples, assign AA=U
                ancestral_allele = 'U'
            elif len(alleles) == 0:
                # If no alleles could be determined from the outgroup samples, assign AA=U
                ancestral_allele = 'U'
            else:
                # Otherwise, the ancestral allele is the one found in the outgroup samples
                ancestral_allele = alleles[0]

            # Update the AA info
            record.info['AA'] = ancestral_allele
            vcf_out.write(record)

```

Script for full filtering:

* First filter is only INFO filtered (QUAL < 20, MQ < 30) **MQ** 
* Second filter is INFO + 5X genotype missingness 10% filtered **5X-MM1** 
* Third filter is above + requiring ancestral polarizable alleles for Cm and Cp. Ancestral alleles assigned as follows: **AA**
  * Take Cm (n=7) and Cp (n=2) from the VCF, if there is only a single allele found among samples, assign that as AA=allele
  * If multiple alleles are found, assign AA=U. Ignore missing genotypes (requiring at least 1 sample with an allele). 
* Last filter is LD pruning, removing any sites with an r2 > 0.2 in 50 SNP windows **LDr2w50** 

It also creates simon martin's invariant sites geno file, using MM1 filters on invariant sites. 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/4.SNP_Filtering.sh ${i}; done 
CHR=$1

mkdir chromosome_vcfs chromosome_vcfs/stats

raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged
chr_map=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/momi2/chr_map.txt

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#Subset samples
bcftools view --threads 10 --samples-file ~/merondun/cuculus_migration/relatedness/Sample_List_Unrelated.list --force-samples -Ou ${raw_vcfs}/${CHR}.vcf.gz | \
        #change chr_1 to 1 
        bcftools annotate --threads 10 --rename-chrs ${chr_map} -Ou | \
        bcftools norm -d both -Oz -o chromosome_vcfs/${CHR}_raw.vcf.gz 
bcftools index --threads 10 chromosome_vcfs/${CHR}_raw.vcf.gz

#filter on DP, retain only sites within mean coverage + 2*sd or mean covearge - 2*sd 
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' chromosome_vcfs/${CHR}_raw.vcf.gz > chromosome_vcfs/${CHR}_raw.dp_stats.txt
mean=$(cat chromosome_vcfs/${CHR}_raw.dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
sd=$(cat chromosome_vcfs/${CHR}_raw.dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
low=$(echo "$mean - 2*$sd" | bc)
high=$(echo "$mean + 2*$sd" | bc)

rm chromosome_vcfs/${CHR}_raw.dp_stats.txt

#subset SNPs, filter for MQ ONLY, fix chromosome name (stripping chr_), so new chromosome is simply 1 
bcftools view --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 --threads 10 -Oz -o chromosome_vcfs/${CHR}_snp.MQ.vcf.gz -i "QUAL > 20 & MQ > 30 & INFO/DP > ${low} & INFO/DP < ${high}" chromosome_vcfs/${CHR}_raw.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ.vcf.gz

#filter for >= 5x genotypes, apply 10% missingness filter
MINDP=5
bcftools +setGT -Ou chromosome_vcfs/${CHR}_snp.MQ.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #update AC fields after changing genotypes to missing 
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -i 'F_MISSING < 0.1' -Oz -o chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz

# Polarize alleles based on C. poliocephalus / C. micropterus shared derived alleles
python ~/merondun/gen_gen/general/Polarize_VCF_Add_AA.py ~/merondun/cuculus_migration/snp_calling/Outgroups.list chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A1

# Switch reference and ancestral allele, ensure --recalc is added to switch AC fields
java -jar ~/modules/jvarkit/dist/vcffilterjdk.jar -f ~/modules/jvarkit/dist/script.js --recalc chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A1.vcf.gz --output chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz

# Remove any alleles which aren't polarizable, and remove genotype fields which weren't switched (e.g. AD) because of ref allele switch
bcftools view --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 --max-af 0.999 --threads 10 -i 'INFO/AA!="U"' chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz -Ob -o - | \
    bcftools annotate -x ^FORMAT/GT,FORMAT/DP,FORMAT/GQ -o chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz
rm chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A1.vcf.gz* chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-A2.vcf.gz*

# LD Prune
bedtools intersect -header -a chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz -b $neutral | \
        bcftools +prune -m 0.2 --window 50 -Oz -o chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz
bcftools index --threads 10  chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz

# And phase with beagle 
java -Xmx40g -jar ~/modules/beagle.28Jun21.220.jar gt=chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz out=chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50-phased nthreads=10 impute=true
bcftools index --threads 10 chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50-phased.vcf.gz

# Grab invariant sites, filter for MQ, DP (same thresholds), and missingness
echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
bcftools view --max-ac 0 -i "MQ > 30 & INFO/DP > ${low} & INFO/DP < ${high} & F_MISSING < 0.1" --threads 10 -Ou chromosome_vcfs/${CHR}_raw.vcf.gz -Oz -o chromosome_vcfs/${CHR}_1N.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}_1N.vcf.gz

# re-merge the invariant and filtered SNPs
bcftools concat --threads 10 -Ob chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz chromosome_vcfs/${CHR}_1N.vcf.gz | \
        bcftools sort -Oz -o chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz
bcftools index --threads 10 chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz

#SNP Counts in tidy format
raw=$(bcftools index -n chromosome_vcfs/${CHR}_raw.vcf.gz)
mq=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ.vcf.gz)
mm1=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ-5X-MM1.vcf.gz)
aa=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz)
ld=$(bcftools index -n chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA-LDr2w50.vcf.gz)
total=$(bcftools index -n chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz)

echo -e "${CHR}\tRaw\t${raw}" > chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tMQ20\t${mq}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\t5X_MM1\t${mm1}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tPolarized\t${aa}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tLD_R2_W50\t${ld}" >> chromosome_vcfs/stats/${CHR}.snp.counts
echo -e "${CHR}\tAllSites\t${total}" >> chromosome_vcfs/stats/${CHR}.snp.counts

#create simon's divergence input file
echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 2 --skipIndels -i chromosome_vcfs/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz | \
        bgzip > chromosome_vcfs/${CHR}.geno.gz

rm chromosome_vcfs/${CHR}_1N.vcf.gz*
```

## Subset N=40

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/7.Subsample_Demography.sh ${i}; done 
CHR=$1

#filtered vcfs directory, polarized alleles
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs
outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs

mkdir $outvcfs

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

#Subset samples
bcftools view --threads 10 --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list -Ov ${vcfs}/${CHR}.MQ-5X-MM1-AA.AllSites.vcf.gz | \
        #RETAIN neutral sites 
        bedtools intersect -header -a - -b $neutral | \
        #EXCLUDE repeats
        bedtools subtract -header -a - -b $repeats | \
        #re-apply 10% missing data threshold 
        bcftools view --threads 10 -i 'F_MISSING < 0.1' -Oz -o ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 10 ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz

#Phase VCF with beagle 
java -jar -Xmx160g ~/modules/beagle.28Jun21.220.jar gt=${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz out=${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp nthreads=8 window=40 overlap=2 impute=true
bcftools index --threads 10 ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz

#add the INFO DP,MQ and FMT/DP annotations back onto this VCF, from the pre-phased VCF 
bcftools annotate --threads 10 -a ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.vcf.gz -c INFO/DP,INFO/MQ,FMT/DP ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz -Oz -o ${vcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.vcf.gz
bcftools index --threads 10 ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.vcf.gz

rm ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz ${outvcfs}/${CHR}.MQ-5X-MM1-AA.AllSites-Neutral-NonRepetitive.Phased.tmp.vcf.gz.csi 
```



## Merge Autosomes 

I will also prune for LD and do a PCA / ADMIXTURE on sample subsets:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# submit as:  for GROUP in $(cat GROUPS.list); do sbatch -J MERGE_${i} ~/merondun/cuculus_migration/admixture/1.Subset_Subgroups.sh ${GROUP}; done 
# mamba activate snps

GROUP=$1

mkdir autosomal_files 

#intergenic and 4-fold degenerate sites 
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

# Merge the VCFs into a single autosomal SNP set
bcftools concat --threads 20 -Ou chromosome_vcfs/*AA.vcf.gz | \
    #subset samples 
    bcftools view --threads 20 -Ou --force-samples --samples-file ~/merondun/cuculus_migration/relatedness/${GROUP}_Unrelated.list | \
    #ensure only variant SNPs are left 
    bcftools view --threads 20 -Ov --min-alleles 2 --max-alleles 2 --min-af 0.05 --max-af 0.999 -i 'F_MISSING < 0.1' | \
    #intersect with neutral sites 
    bedtools intersect -header -a - -b $neutral | \
    #prune LD 
    bcftools +prune -m 0.2 --window 50 -Oz -o autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.vcf.gz
bcftools index --threads 20 autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.vcf.gz

#Create plink bed files for admixture + perform pca 
~/modules/plink2 --threads 20 --vcf autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.vcf.gz --chr-set 29 --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --recode vcf bgz --pca --out autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50
```

Create PCA and run ADMIXTURE on the full dataset:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# submit as: for GROUP in $(cat GROUPS.list); do for K in {2..10}; do sbatch -J ADMIX_${GROUP}_${K} ~/merondun/cuculus_migration/admixture/2.ADMIXTURE.sh ${GROUP} ${K}; done ; done

GROUP=$1
K=$2

mkdir -p admixture/unrelated_n261/${GROUP}
cd admixture/unrelated_n261/${GROUP}

#Run Admixture
admixture -j7 --cv=5 ../../autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50.bed ${K} > ${GROUP}.MQ-5X-MM1-AA-LDr2w50.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../../autosomal_files/${GROUP}.MQ-5X-MM1-AA-LDr2w50 -fname ${GROUP}.MQ-5X-MM1-AA-LDr2w50.${K}.P -qname ${GROUP}.MQ-5X-MM1-AA-LDr2w50.${K}.Q -P 10 -o eval_${GROUP}.MQ-5X-MM1-AA-LDr2w50_${K}
```

# Analyses

## Ancestral Reconstruction

Grab 5K gene trees with 6670 OTUs using the 'Sequenced Species' dataset of Hackett from [BirdTree](https://birdtree.org/subsets/)

Includes 126 of the Cuculiformes species:

```bash
Stage 2 MayrParSho Hackett: 10K trees with 6670 OTUs [10k sampled]
Job ID: tree-pruner-3972b4dc-d4e3-42db-8473-30777033c335
```

Estimate species tree from the hackett trees:

```bash
# First estimate species tree
java -jar ~/modules/ASTRAL/astral.5.7.7.jar -i stage2_mayrparshohackett_10k.tre -o stage2_mayrparshoahackett_10k_species_tree.tre

# Estimate concordance factors using the species tree against the 10K hackett trees
iqtree -t stage2_mayrparshohackett_10k_species_tree.tre --gcf stage2_mayrparshohackett_10k.tre --prefix stage2_mayrparshohackett_10k_concord
```

| Route | Nonbreeding | Description                                      |
| ----- | ----------- | ------------------------------------------------ |
| S     | I           | Sedentary within Indo-Malay / Australia          |
| S     | SA          | Sedentary within the Americas                    |
| S     | A           | Sedentary within Africa                          |
| II    | I           | Migration within Indo-Malay / Australia          |
| SS    | SA          | Migration within South America                   |
| AA    | A           | Migration within Africa                          |
| APw   | A           | Africa - Palearctic migration, western route     |
| APe   | A           | Africa - Palearctic migration, eastern route     |
| APwe  | A           | Africa - Palearctic migration, polymorphic route |
| NS    | SA          | North - South America migration                  |

Using that concordance tree, estimate ancestral state changes:

```bash
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
pdf('20250115_ReconstructionER-Cuculiformes_Residency-Node80.pdf',height=7,width=7)
cols<-setNames(viridis(3)[1:length(unique(phenotypes))],sort(unique(phenotypes)))
plotTree(t2,ftype="off")

#Node points
node_labs <- ifelse(round(as.numeric(t2$node.label)) >= 80, NA, '*')
nodelabels(node_labs,node=1:t2$Nnode+Ntip(t2),cex=0.7,bg=NA,adj = c(0.5,-.5))
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
pdf('20250115_ReconstructionER-Cuculiformes_Route-Node80-Tall.pdf',height=7,width=3.5)
cols<-setNames(viridis(7,option='turbo')[1:length(unique(phenotypes2))],sort(unique(phenotypes2)))
cols <- setNames(c('#327739','#b5d887','#5db465','darkorchid2','#b54325','white','grey90','grey50')[1:length(unique(phenotypes2))],sort(unique(phenotypes2)))
t2_clado <- compute.brlen(t2, method = "equal")  # Assign equal branch lengths
plotTree(t2_clado,ftype="off")

tip_colors2 <- cols[phenotypes2] 
tiplabels(pch=22, bg=tip_colors2, col="black", cex=0.5,adj=c(0,0))

node_labs <- ifelse(round(as.numeric(t2$node.label)) >= 80, NA, '*')
nodelabels(node_labs, node=1:t2$Nnode+Ntip(t2), cex=1, bg=NA, adj=c(0.5, 0), frame="none")
nodelabels(node=1:t2$Nnode+Ntip(t2),
           pie=fitER2$lik.anc,piecol=cols,cex=0.8)
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

```







## ADMIXTURE + Tesselation Figures

Resubsetting samples:

```bash
GROUP=Canorus
GROUP=Optatus

~/modules/plink2 \
    --bfile CanorusOptatus-n261.MQ-5X-MM1-AA-LDr2w50 \
    --chr-set 29 \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --keep ~/merondun/cuculus_migration/admixture/${GROUP}.pid \
    --geno 0.10 \
    --maf 0.05 \
    --make-bed \
    --recode vcf bgz \
    --pca \
    --out ${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50
```

Ensure sample names don't have any weird plink changes:

```bash
bcftools reheader --samples N261.renamed --threads 10 --output CanorusOptatus-n261.MQ-5X-MM1-AA-LDr2w50.renamed.vcf.gz CanorusOptatus-n261.MQ-5X-MM1-AA-LDr2w50.vcf.gz
bcftools query -l CanorusOptatus-n261.MQ-5X-MM1-AA-LDr2w50.vcf.gz | sed 's/_M_.*/_M/g' | sed 's/_F_.*/_F/g' > N261.renamed
```

Run admixture:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# mamba activate snps
# submit as: for GROUP in $(cat GROUPS.list); do for K in {2..10}; do sbatch -J ADMIX_${GROUP}_${K} ~/merondun/cuculus_migration/admixture/2.ADMIXTURE.sh ${GROUP} ${K}; done ; done

GROUP=$1
K=$2

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/unrelated_n261/${GROUP}
bedfiles=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/autosomal_files

mkdir -p ${wd}
cd ${wd}

#Run Admixture
admixture -j7 --cv=5 ${bedfiles}/${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.bed ${K} > ${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ${bedfiles}/${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50 -fname ${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.${K}.P -qname ${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50.${K}.Q -P 10 -o eval_${GROUP}-n261.MQ-5X-MM1-AA-LDr2w50_${K}

```

The follow R script will make these plots:

* evalAdmix plots

* structure plots

* tesselation plots for K = 2

* tesselation plots with pie charts for K  = 2 - 5

* blank map plot with the study extent

* map with the sampling of the population genetic samples

  ````
  ### Plot ADMIXTURE
  setwd('~/merondun/cuculus_migration/admixture/full_qs/')
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
  library(RColorBrewer)
  library(ggnewscale)
  library(scales)
  
  prefixes = gsub('.2.Q','',list.files('.',pattern='.*.2.Q'))
  prefixes = c("Canorus-n261.MQ-5X-MM1-AA-LDr2w50","Optatus-n261.MQ-5X-MM1-AA-LDr2w50")
  qdir = '.' #directory with Q files
  counter = 0 
  cols = brewer.pal(8,'Paired')[c(3,4,7,8)]
  
  for (admix_run in prefixes) { 
    counter = counter + 1
    cat('Working on run: ',admix_run,'\n')
    
    #Set if else parameters for colors, species, and shape based on the prefix string 
    if(grepl('Canorus',admix_run)){ 
      colz=cols[c(2,1)]; sp='CC'; spshape=21; filt='Cuculus canorus'
    } else {
      colz=cols[c(3,4)]; sp='CO'; spshape=24; filt = 'Cuculus optatus'
    }
    
    prefix = admix_run 
    admix = melt_admixture(prefix = prefix, qdir = qdir)
    
    #read in metadata
    md = read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')
    
    # Change optatus east to just gray
    md <- md %>% mutate(Population = ifelse(Population == 'COE',NA,Population),
                        PopColor = ifelse(Population=='COE','grey60',PopColor))
    
    # Assign groups based on the n=40 demographic inference subset
    n40 = read_tsv('~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.pop',col_names = F)
    names(n40) = c('ID','Group')
    md = md %>% left_join(.,n40) %>% mutate(GroupBinary = ifelse(is.na(Group),'Other','Target'),
                                            Group = ifelse(is.na(Group),'Other',Group))
    
    # Only retain unrelated samples, no toepads
    md_sub = md %>% filter(SpeciesShort == sp)
    
    #md = read_tsv('~/merondun/cuculus_migration/ebird/Metadata_Breeding_2024APR24.txt')
    admixmd = left_join(md_sub,admix)
    
    # Reorder individuals baseed on longitude
    admixmd = admixmd %>% mutate(ID = fct_reorder(ID,Longitude))
    
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
      scale_x_continuous(breaks = seq(min(evaldat$K), max(evaldat$K), by = 1)) +
      coord_flip()
    ep
    assign(paste0('e',1),ep)
    
    ggsave(paste0('~/merondun/cuculus_migration/figures/20250115_evalAdmix_',sp,'.pdf'),e1,height=5,width=6,dpi=600)
    
    #Now plot ADMIXTURE 
    adplot =
      admixmd %>% filter(Specified_K == 2 ) %>%  #specify the levels you want 
      mutate(Specified_K = paste0('K',Specified_K)) %>% 
      ggplot(aes(x = factor(ID), y = Q, fill = factor(K), col=GroupBinary)) +
      geom_col(size = 0.1) +
      facet_grid(Specified_K~Species, scales = "free", space = "free") +
      theme_minimal(base_size=6) + labs(x = "",y = "") +
      scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
      scale_fill_manual(values=colz)+
      scale_color_manual(values=c('gray','black'))+
      theme(
        panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x = element_blank(),
        axis.text.x=element_text(angle=90,size=5),
        axis.text.y = element_text(size=3),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(size=6)
      )
    adplot
    assign(paste0('p',counter),adplot)
    
    ggsave(paste0('~/merondun/cuculus_migration/figures/20250115_Admixture_',sp,'.pdf'),adplot,height=2.5,width=6,dpi=600)
    
    #save the K1/K2 Q values:
    qval = admixmd %>% filter(Specified_K == 2 & K == 'K1') %>% select(ID,K1 = Q)
    write.table(qval,paste0('~/merondun/cuculus_migration/admixture/Assigned_K2_Qvalues_',prefix,'.txt'),quote=F,sep='\t',row.names=F)
    
    ### Tesselation
    # Import birdlife shapefiles, breeding == 2, presence ==1 means extant
    bg <-  st_read('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/spatial/birdlife_breeding_distributions/SppDataRequest.shp')
    filtered_data <- bg[bg$PRESENCE == 1 & bg$SEASONAL == 2 & bg$SCI_NAME == filt, ]
    
    #which K value to plot 
    show_k = 2
    
    # Import Q, and then merge with the fam file to get the IDs to retain
    show_q = read.table(paste0(qdir,'/',admix_run,'.',show_k,'.Q')) #read in the specific file
    fam = read_tsv(paste0(prefix,'.fam'),col_names = F) %>% select(X2) %>% dplyr::rename(ID = X2)
    all_q = cbind(show_q,fam)
    
    # This will only retain the individuals from the admixture lot 
    retain_input = left_join(admixmd %>% select(ID,Group,Latitude,Longitude) %>% unique,all_q)
    show_q_mat = as.matrix(retain_input %>% select(V1,V2)) #convert it to a matrix
    class(show_q_mat) = c('tess3Q','matrix','array') #make sure tess3r thinks that it's actually a tess object
    coords = retain_input %>% select(Longitude,Latitude) #convert lat and long
    coords_mat = as.matrix(coords) #convert coordinates to matrix
    
    #plot using ggplot
    map.polygon <- getMap(resolution = "high")
    # Jitter the lat/long 
    sitesp = st_as_sf(retain_input %>% select(ID,Group,Longitude,Latitude) %>% distinct %>% 
                        mutate(loj = jitter(Longitude,amount=0.5),laj = jitter(Latitude,amount=0.5)),remove = F, coords = c("loj", "laj"), crs = 4326, agr = "constant") 
    sitesp = sitesp %>% mutate(Kept = ifelse(Group == 'Other','Stronghold','Other'))
    pl = ggtess3Q(show_q_mat, coords_mat, map.polygon = filtered_data,col.palette = colz)
    k2p1 = pl +
      geom_path(data = map.polygon, aes(x = long, y = lat, group = group),col='white',lwd=0.2) +
      coord_sf(xlim = c(-20,185), 
               ylim = c(-30,75), expand = FALSE)+  
      new_scale_fill()+
      new_scale_color()+
      geom_point(data = sitesp, aes(x = loj, y = laj,fill=Group,alpha=Kept),shape=spshape, size = 2,col='black') +
      scale_fill_manual(values=md$PopColor,breaks=md$Population)+
      scale_alpha_manual(values=c(0.8,0.2))+
      xlab("Longitude") + ylab("Latitude") + 
      theme_classic(base_size=8)+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
      theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
    k2p1
    
    ggsave(paste0('~/merondun/cuculus_migration/figures/20250115_Tesselation_',sp,'.pdf'),k2p1,height=4,width=7,dpi=300)
    
    #plot all K...
    maximum_k = 5
    for (kval in seq(2,maximum_k,1)) {
      
      cat('Working on K = ',kval,'\n')
      # Import Q, and then merge with the fam file to get the IDs to retain
      show_q = read.table(paste0(qdir,'/',admix_run,'.',kval,'.Q')) #read in the specific file
      fam = read_tsv(paste0(prefix,'.fam'),col_names = F) %>% select(X2) %>% dplyr::rename(ID = X2)
      all_q = cbind(show_q,fam)
      
      # This will only retain the individuals from the admixture lot 
      retain_input = left_join(admixmd %>% select(ID,Group,Latitude,Longitude) %>% unique,all_q)
      show_q_mat = as.matrix(retain_input %>% select(matches('V'))) #convert it to a matrix
      class(show_q_mat) = c('tess3Q','matrix','array') #make sure tess3r thinks that it's actually a tess object
      coords = retain_input %>% select(Longitude,Latitude) #convert lat and long
      coords_mat = as.matrix(coords) #convert coordinates to matrix
      
      # First, grab the individuals and calculate the mean Q values within each cluster. Cluster will be geographic reigon, also calculate mean lat/long for plotting
      kept = admixmd %>% select(ID,Specified_K,K,Q,Latitude,Longitude,Group = GeographicGroup)
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
        ggtess3Q(show_q_mat, coords_mat, map.polygon = filtered_data,col.palette = viridis(kval)) + 
        geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
        geom_sf(data = sites, aes(geometry = geometry), size = 0.1, alpha = 0.1, pch=26) +
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
    
    png(paste0('~/merondun/cuculus_migration/figures/20250115_Tesselation_',admix_run,'_AllK.png'),res=600,units='in',height=7,width=10)
    print(ggarrange(t2,t3,t4,t5,ncol=2,nrow=2))
    dev.off()
    
  }
  
  
  ###### Blank plot ######
  #Create a small psuedo region for plotting the tesselation, makes a small file size 
  coords <- matrix(c(50, 60, 51, 60, 51, 61, 50, 61, 50, 60), ncol = 2, byrow = TRUE)
  polygon <- st_polygon(list(coords))
  sf_df <- st_sf(geometry = st_sfc(polygon))
  
  #plot using ggplot
  map.polygon <- getMap(resolution = "high")
  b = ggtess3Q(show_q_mat, coords_mat, map.polygon = sf_df,col.palette = rev(cols))
  blank = b +
    geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='grey90',lwd=0.2) +
    coord_sf(xlim = c(-20,185), 
             ylim = c(-30,75), expand = FALSE)+  
    new_scale_fill()+
    new_scale_color()+
    geom_point(data = sitesp, aes(x = loj, y = laj,fill=Group,alpha=Kept),shape=26,size = 2,col='black') +
    xlab("Longitude") + ylab("Latitude") + 
    theme_classic(base_size=8)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
    theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
  blank
  
  ggsave(paste0('~/merondun/cuculus_migration/figures/20241023_Tesselation_Blank.pdf'),blank,height=4,width=7,dpi=300)
  
  #### Add individual points ####
  sitesp = st_as_sf(md %>% filter(Analysis_ADMIXTURE_All == 1) %>% select(ID,Group,Longitude,Latitude,SpeciesShort) %>% distinct %>% mutate(loj = jitter(Longitude,amount=1),laj = jitter(Latitude,amount=1)),remove = F, coords = c("loj", "laj"), crs = 4326, agr = "constant") 
  sitesp = sitesp %>% mutate(Kept = ifelse(Group == 'Other','Stronghold','Other'))
  pos = b +
    geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='grey90',lwd=0.2) +
    coord_sf(xlim = c(-20,185), 
             ylim = c(-30,75), expand = FALSE)+  
    new_scale_fill()+
    new_scale_color()+
    geom_point(data = sitesp, aes(x = loj, y = laj,fill=Group,alpha=Kept,shape=SpeciesShort), size = 2,col='black') +
    scale_fill_manual(breaks=c('CCW','CCE','COW','COE','Other'),values=c("#1F78B4","#A6CEE3","#E31A1C","#FB9A99","white" ))+
    xlab("Longitude") + ylab("Latitude") + 
    scale_shape_manual(values=c(21,24))+
    guides(fill=guide_legend(nrow=2,override.aes=list(shape=22)))+
    scale_alpha_manual(values=c(0.8,0.4))+
    theme_classic(base_size=8)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
    theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
  pos
  
  pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/PopGenLocations_Map_Tesselation_2024MAY14.pdf',height=4,width=7)
  pos
  dev.off()
  ````
  
  

![image-20250122103827205](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20250122103827205.png)

Structure plots for canorus and optatus, run separately

![image-20250122104134393](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20250122104134393.png)

Inferred tesselations of ancestry projected across geographic space for K=2, separately for canorus and optatus, and clipped to birdlife international's breeding extant range.

![image-20250122104905746](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20250122104905746.png)

Tesselations for K2 - K5 with pie charts showing proportions of ancestry within each region

## PCA

```R
setwd('~/merondun/cuculus_migration/admixture/')
library(tidyverse)
library(RColorBrewer)

# Reorder individuals baseed on longitude
md <- read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')
md <- md %>% mutate(Population = ifelse(Population == 'COE',NA,Population),
                    PopColor = ifelse(Population=='COE','grey60',PopColor))

groups = c('CanorusOptatus','Canorus','Optatus')

for (group in groups) {
  
  #read PCA data
  vec = read.table(paste0(group,'-n261.MQ-5X-MM1-AA-LDr2w50.eigenvec'),header=TRUE,comment.char = '')
  vec = vec %>% dplyr::rename(ID = X.FID) %>% select(-IID)
  val = read.table(paste0(group,'-n261.MQ-5X-MM1-AA-LDr2w50.eigenval'),header=FALSE,comment.char = '')
  val = val %>% mutate(VE = paste0(round(V1/sum(V1),2)*100,'%'))
  
  admixmd = left_join(vec,md) %>% mutate(ID = fct_reorder(ID,Longitude)) %>% mutate(Population = ifelse(is.na(Population),'Unassigned',Population))
  
  pca_plot <- admixmd %>% 
    ggplot() +
    geom_point(aes(x=PC1,y=-PC2,fill = Population,shape=Species_Latin,alpha=as.factor(Analysis_Demography)), 
               size = 2, color = 'grey20') +
    theme_bw(base_size=8) +
    scale_fill_manual(values=md$PopColor,breaks=md$Population)+
    scale_shape_manual(values=md$Shape,breaks=md$Species_Latin)+
    scale_alpha_manual(values=c(0.25,0.9))+
    xlab(paste0('PC1 (',val[1,2],')'))+
    ylab(paste0('PC2 (',val[2,2],')'))+
    guides(fill = guide_legend(override.aes = list(shape = 21)))+
    coord_flip()
  pca_plot
  
  ggsave(paste0('../figures/20250115_PCA_',group,'.pdf'),pca_plot,height=2,width=4,dpi=300)
         
}

```

![image-20241023160342428](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20241023160342428.png)

## FST ~ Geographic Distance

```R
#### Calculate pairwise geographic distance between distance groups 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/isolation_by_distance')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(sf)
library(ggspatial)
library(RColorBrewer)

#Read in metadata
md <- read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')

#groups for fst: distance ~  autos
analyze_these = md %>% 
  filter(Analysis_ADMIXTURE_All == 1) %>% 
  select(GeographicGroup) %>% 
  group_by(GeographicGroup) %>% 
  mutate(DistanceHaps = n()) %>% 
  unique %>% 
  filter(DistanceHaps >= 3) %>% pull(GeographicGroup)

#calculate geographic distance between the groups
dists <- md %>% filter(Analysis_ADMIXTURE_All == 1 & GeographicGroup %in% analyze_these) %>% group_by(GeographicGroup) %>% summarize(meanLat = mean(Latitude),meanLong = mean(Longitude))

dists <- dists %>% mutate(GCol = rep_len(c(viridis(9, option = 'turbo'),rev(viridis(9, option = 'turbo'))), 21), 
                 GShape = rep_len(c(21, 24, 25, 22), 21))
dists %>% select(GCol,GShape) %>% unique %>% nrow

#plot 
world = map_data("world")
sites = st_as_sf(dists, coords = c("meanLong", "meanLat"), 
                 crs = 4326, agr = "constant") 
fst_compars = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  geom_sf(data = sites, 
          aes(fill=GeographicGroup,shape=GeographicGroup),
          size=4,alpha=0.9,show.legend = T,stroke=0.5) +
  scale_shape_manual(values=dists$GShape,breaks=dists$GeographicGroup)+
  scale_fill_manual(values=dists$GCol,breaks=dists$GeographicGroup)+
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(dists$meanLong)-5, max(dists$meanLong)+5), 
           ylim = c(min(dists$meanLat)-5, max(dists$meanLat)+5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)
fst_compars

ggsave('../figures/20250121_FST_IBD_Populations.pdf',
       fst_compars,height=4,width=9,dpi=600)

# initialize an empty data frame to store pairwise distances
pairwise_distances <- data.frame()
sub_data = dists

# calculate pairwise distances for the current 'W' GeographicGroup
for(i in 1:(nrow(sub_data) - 1)) {
  for(j in (i + 1):nrow(sub_data)) {
    
    point1 <- c(sub_data$meanLong[i], sub_data$meanLat[i])
    point2 <- c(sub_data$meanLong[j], sub_data$meanLat[j])
    
    distance_km <- distHaversine(point1, point2) / 1000  # convert to km
    
    # append the result to the pairwise_distances data frame
    pairwise_distances <- rbind(pairwise_distances, 
                                data.frame(P1 = sub_data$GeographicGroup[i], 
                                           P2 = sub_data$GeographicGroup[j],
                                           Distance_km = distance_km))
  }
}

pairwise_distances = pairwise_distances %>% 
  mutate(Group = paste0(P1,'__',P2)) %>% 
  separate(P1,into=c('S1','G1'),remove=F) %>% 
  separate(P2,into=c('S2','G2'),remove=F) %>% 
  filter(S1 == S2) %>% select(P1,P2,Distance_km,Group)


# show the calculated pairwise distances
write_tsv(pairwise_distances,file='Pairwise_GeographicDistance_Km.txt')

#write out populations
for (pop in unique(analyze_these)) {
  su = md %>% filter(Analysis_ADMIXTURE_All == 1 & GeographicGroup == pop)
  write.table(su$ID,file=paste0('populations/',pop,'.pop'),quote=F,sep='\t',row.names=F,col.names=F)
}

```

![image-20250121112454635](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20250121112454635.png)

Comparisons:

```bash
awk '{print $4}' Pairwise_GeographicDistance_Km.txt | sed '1d' > PairwiseComparisons.list
cat PairwiseComparisons.list 
CC_5__CC_9
CC_6__CC_7
CC_6__CC_9
CC_7__CC_9
CO_10__CO_2
CO_10__CO_4
CO_10__CO_7
```

Calculate FST:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

# mamba activate snps
# for i in $(cat PairwiseComparisons.list); do sbatch -J FST_${i} 2.Pairwise_Distance_FST.sh ${i}; done 
GROUP=$1

VCF=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/autosomal_files/CanorusOptatus-n261.MQ-5X-MM1-AA-LDr2w50.vcf.gz
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/isolation_by_distance

for CHR in $(cat Chromosomes.list); do

echo "WORKING ON CHR: ${GROUP} and ${CHR}"

p1=$(echo ${GROUP} | sed 's/__.*//g')
p2=$(echo ${GROUP} | sed 's/.*__//g')

mkdir -p work out

#calculate fst
~/modules/vcftools/bin/vcftools --gzvcf ${VCF} --fst-window-size 100000 --fst-window-step 100000 --weir-fst-pop populations/${p1}.pop --weir-fst-pop populations/${p2}.pop --out work/${CHR}_${GROUP}
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3,$4,$5, g}' work/${CHR}_${GROUP}.windowed.weir.fst > out/${CHR}_${GROUP}.fst

fi

done
```

Merge output: `mergem *fst > ../20250122_Pairwise_FST_DistanceGroups_Autosomes.txt`

Plot final correlations:

```R
#### Plot Distance ~ FST
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/isolation_by_distance')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(broom)
library(ggpubr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(meRo)
library(RColorBrewer)

#Read in metadata
md = read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')

d = read_tsv('20250122_Pairwise_GeographicDistance_Km.txt')
f = read_tsv('20250122_Pairwise_FST_DistanceGroups_Autosomes.txt')

#prep and merge frames 
names(f) = c('chr','start','end','snps','FST','Group')
f2 = f %>% mutate(FST = pmax(0, pmin(1, FST))) %>% select(!c(start,end,snps))

#merge with geographic distance 
df = left_join(f2,d) 

dfs = df %>% group_by(Group,P1,P2,Distance_km) %>% 
  sum_stats(FST)
dfs %>% ggplot(aes(x=log(Distance_km),y=log(mean)))+
  geom_point()+
  geom_smooth()+
  theme_bw()

dfs <- dfs %>% mutate(Species = ifelse(grepl('CO',Group),'C. optatus','C. canorus'))

# Spearman's rho and cor test within each AvZ level and gather coefficients
cors = dfs %>% group_by(Species) %>%
  summarize(
    rho = cor(mean, Distance_km,method='spearman'),
    p_value = cor.test(mean, Distance_km,method='spearman')$p.value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'),
         signif = case_when(
           padj < 0.05 ~ "*",
           TRUE ~ " (n.s.)"),
         label = paste0('rs: ',round(rho, 2), signif))
cors
# Species      rho p_value  padj signif    label          
# <chr>      <dbl>   <dbl> <dbl> <chr>     <chr>          
#   1 C. canorus 0.772   0     0     "*"       rs: 0.77*      
#   2 C. optatus 0.221   0.427 0.853 " (n.s.)" rs: 0.22 (n.s.)

cols = c('darkorchid2','orange')
dfs$Species = factor(dfs$Species,levels=c('C. canorus','C. optatus'))

#lmm, account for P1 and P2 
model_summaries <- dfs %>%
  group_by(Species) %>%
  do({
    model <- lmer(log(pmax(mean, 0.005)) ~ log(pmax(Distance_km, 0.005)) + (1|P1) + (1|P2), data = .)
    summary_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  }) %>%
  dplyr::bind_rows()

#output fixed effect of distance, bonferroni correction 
model_summaries %>% filter(grepl('Distance',term)) %>% select(-effect,-term) %>% ungroup %>% 
  mutate(padj = p.adjust(p.value,method='bonferroni'),
         signif = ifelse(padj < 0.05,'*','n.s.'))

# Species    estimate std.error statistic    df  p.value conf.low conf.high     padj signif
# <fct>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr> 
#   1 C. canorus    0.743    0.0465     16.0  98.8  3.82e-29    0.651     0.836 7.64e-29 *     
#   2 C. optatus    0.446    0.0691      6.45  5.43 9.77e- 4    0.272     0.619 1.95e- 3 *  


# Create the plot
pp1 = ggplot(dfs, aes(x = log(pmax(Distance_km,0.005)), y = log(pmax(mean,0.005)), color = Species, shape = Species)) +
  geom_point(size=1) + #0.25 for main plot 
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +  #0.5 for main plot 
  geom_text(data = cors[1,], aes(x = -Inf, y = Inf, label = paste0('C. canorus rho: ',round(rho,3))), size = 3,vjust = 2, hjust = -.25) +
  geom_text(data = cors[2,], aes(x = -Inf, y = Inf, label = paste0('C. optatus rho: ',round(rho,3))), size = 3,vjust = 4, hjust = -.25) +
  labs(title = "Relationship between FST and Distance_km",
       x = "log(Geographic Distance (km))",
       y = "log(FST)") +
  scale_color_manual(values=cols)+
  theme_bw(base_size=6) + theme(legend.position = 'none')
pp1


ggsave('~/merondun/cuculus_migration/figures/20250115_FSTvGeoDistanceSpecies.pdf',pp1,height=3,width=3.5,dpi=300)


#### Fst / (1 - FST) and log(dist)

cors_rousset = dfs %>% group_by(Species) %>%
  summarize(
    rho = cor(mean, Distance_km,method='spearman'),
    p_value = cor.test( (mean / (1 - mean)), Distance_km,method='spearman')$p.value) %>% 
  mutate(label = paste0('rs: ',round(rho, 2)))
cors_rousset

#lmm, account for P1 and P2 
model_summaries_rousset <- dfs %>%
  group_by(Species) %>%
  do({
    model <- lmer( (mean / (1 - mean)) ~ log(pmax(Distance_km, 0.005)) + (1|P1) + (1|P2), data = .)
    summary_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  }) %>%
  dplyr::bind_rows()

#output fixed effect of distance, bonferroni correction 
model_summaries_rousset %>% filter(grepl('Distance',term)) %>% select(-effect,-term) %>% ungroup %>% 
  mutate(padj = p.adjust(p.value,method='bonferroni'),
         signif = ifelse(padj < 0.05,'*','n.s.'))

# Species    estimate std.error statistic    df  p.value conf.low conf.high     padj signif
# <fct>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr> 
#   1 C. canorus  0.0310   0.00220      14.1  89.9  1.85e-24  0.0266    0.0354  3.70e-24 *     
#   2 C. optatus  0.00737  0.000740      9.96  5.14 1.49e- 4  0.00548   0.00926 2.98e- 4 *   


# Create the plot
pp2 = ggplot(dfs, aes(x = log(pmax(Distance_km,0.005)), y = (mean / (1 - mean)), color = Species, shape = Species)) +
  geom_point(size=1) + #0.25 for main plot 
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +  #0.5 for main plot 
  geom_text(data = cors[1,], aes(x = -Inf, y = Inf, label = paste0('C. canorus rho: ',round(rho,3))), size = 3,vjust = 2, hjust = -.25) +
  geom_text(data = cors[2,], aes(x = -Inf, y = Inf, label = paste0('C. optatus rho: ',round(rho,3))), size = 3,vjust = 4, hjust = -.25) +
  labs(title = "Relationship between FST and Distance_km",
       x = "log(Geographic Distance (km))",
       y = "FST / (1 - FST)") +
  scale_color_manual(values=cols)+
  theme_bw(base_size=6) + theme(legend.position = 'none')
pp2


ggsave('~/merondun/cuculus_migration/figures/20250115_FSTvGeoDistanceRoussetSpecies.pdf',pp1,height=3,width=3.5,dpi=300)


```



## Intersect with Migration Data

```bash
### Plot ADMIXTURE
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/unrelated_all/')
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
library(RColorBrewer)
library(ggnewscale)
library(geosphere)

#### Identify western and eastern strongholds #### 

#we want to assign migratory tracks based on their breeding locations, so first, let's import the K1/K2 ancestry coefficients and find the longitude of the most 100% western and eastern groups as our threshold
# canorus: K1 = WEST, optatus: K1 = EAST
md <- read_tsv('~/merondun/cuculus_migration/Full_Metadata.txt')
#create a bounding box for each species based on the K1 and K2 strongholds
thresholds = md %>% 
  drop_na(K1) %>%
  # Assign samples based on their qvalue and species 
  mutate(Group = ifelse(K1 > 0.99 & SpeciesShort == 'CC','CCW',
                        ifelse(K1 < 0.01 & SpeciesShort == 'CC','CCE',
                               ifelse(K1 > 0.99 & SpeciesShort == 'CO','COE',
                                      ifelse(K1 < 0.01 & SpeciesShort == 'CO','COW','Unassigned'))))) %>% 
  filter(Group != 'Unassigned') %>% 
  group_by(Group) %>% 
  summarize(longmn = min(Longitude),
            longmx = max(Longitude),
            latmn = min(Latitude),
            latmx = max(Latitude))
thresholds

# Plot those regions 
map.polygon <- getMap(resolution = "high")
strong = ggplot() + 
  geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='grey90',lwd=0.2) +
  coord_sf(xlim = c(-20,185), 
           ylim = c(-30,75), expand = FALSE)+  
  new_scale_fill()+
  new_scale_color()+
  geom_rect(data = thresholds, aes(xmin = longmn, xmax = longmx, ymin=latmn, ymax=latmx, fill=Group),alpha=0.5,col='black',lty=2) +
  scale_fill_manual(breaks=c('CCW','CCE','COW','COE','Other'),values=c("#1F78B4","#A6CEE3","#E31A1C","#FB9A99","white" ))+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
png('../../figures/ADMIXTURE_StrongholdsK2_2024MAY15.png',units='in',res=300,height=4,width=7)
strong
dev.off()

#### Add Migration Data ####
m <- read_csv('../../spatial/03-BoD-deaths-removed-no-NDVI.csv')
m %>% select(dataset, individual.taxon.canonical.name,tag.local.identifier) %>% unique %>% count(individual.taxon.canonical.name)

tracks_df = m %>% filter(!grepl('satura',individual.taxon.canonical.name)) %>% select(individual.taxon.canonical.name,timestamp,Route,Wolf,tag.local.identifier,location.long,location.lat)
names(tracks_df) = c('Species','DT','Route','WolfID','ID','Long','Lat')

#assign 'year1', 'year2', 'year3', etc, starting with 'november' as the base
tracks_df2 = tracks_df %>%
  group_by(ID) %>%
  mutate(DT = mdy_hm(DT)) %>% 
  arrange(DT) %>%
  # For each cuckoo, assign a year identifier - the 'year' starts at November as they were rung in the summer, so summer - nov = year 0 
  # This will also add the number of days in between GPS locations
  mutate(
    Year_Num = cumsum(month(DT) == 11 & lag(month(DT), default = 10) != 11),  # Increment each November
    Year_Label = paste0("Year", Year_Num),
    Days_Between = as.numeric(DT - lag(DT),units='days')
  ) %>%
  ungroup() %>%
  # Also add how many years were sampled for each cuckoo 
  # This will also assign the geographic distance (in km) between consecutive GPS positions
  group_by(ID) %>% 
  mutate(Years = n_distinct(Year_Label),
         Prev_Long = lag(Long),
         Prev_Lat = lag(Lat),
         Distance = (distHaversine(cbind(Prev_Long, Prev_Lat), cbind(Long, Lat))/1000)) %>% 
  ungroup %>% select(-Prev_Long,-Prev_Lat) %>% 
  # And finally, add a case condition for if the points are during a typical breeding or non-breeding season 
  mutate(
    Status = case_when(
      month(DT) %in% c(1,12) ~ "Overwintering",
      month(DT) %in% c(6) ~ "Breeding",
      #(month(DT) == 5 & day(DT) >= 15) | (month(DT) == 6 & day(DT) <= 15) ~ "Breeding",
      TRUE ~ NA_character_
    )) %>% 
  # Replace NA in distance / days_between with 0 since they are the initial point
  replace_na(list(Days_Between = 0, Distance = 0)) 

# Plot time on the x axis, with the distance between points on the y axis for all indivdiauls by color to see if there are migration peak times
monthdist = tracks_df2 %>%
  ggplot(aes(x=month(DT),y=Distance,col=ID))+
  geom_point()+
  geom_line()+
  xlab('Month')+ylab('Distance between fixes (km)')+
  theme_bw()+
  theme(legend.position='none')
png('../../figures/Migration_MonthlyDistances_2024MAY15.png',units='in',res=300,height=3,width=6)
monthdist
dev.off()

# # Distance traveled by month
meandist = tracks_df2 %>%
  group_by(month(DT)) %>%
  sum_stats(Distance) %>%
  ggplot(aes(x=as.factor(`month(DT)`),y=mean,ymin=conf_low,ymax=conf_high))+
  geom_point()+ ylab('Sequential Distance (km)')+xlab('Month')+
  geom_errorbar()+
  theme_bw()
png('../../figures/Migration_MeanMonthlyDistances_2024MAY15.png',units='in',res=300,height=3,width=6)
meandist
dev.off()

# Summaries of distance between points and time between points
tracks_df2 %>% 
  select(Species,ID,Days_Between,Distance) %>% 
  pivot_longer(c(Days_Between,Distance)) %>% 
  ggplot(aes(x=value,fill=Species))+
  geom_histogram()+
  facet_wrap(name~.,scales='free')+
  theme_bw()

#here, simply assign mean breeding as mean location over all years in june 
breed_dat <- tracks_df2 %>% 
  #filter(Days_Between < 10 & Distance < 1000) %>% 
  #filter(Status == 'Breeding') %>% 
  #filter(Status == 'Overwintering') %>% 
  drop_na(Status) %>% 
  group_by(Species,ID,Status) %>% 
  # Assign mean lat and long during the breeding season 
  summarize(Longitude = mean(Long), Latitude=mean(Lat)) %>% 
  ungroup

# Function to check if coordinates fall within a bounding box
in_bbox <- function(longitude, latitude, bbox) {
  longitude >= bbox$longmn & longitude <= bbox$longmx & latitude >= bbox$latmn & latitude <= bbox$latmx
}

# Assign breeding locations to groups based on species and bounding boxes
breed_dat <- breed_dat %>%
  mutate(Group = case_when(
    Species == "Cuculus canorus" & Status == "Breeding" & in_bbox(Longitude, Latitude, thresholds[thresholds$Group == "CCE", ]) ~ "CCE",
    Species == "Cuculus canorus" & Status == "Breeding" & in_bbox(Longitude, Latitude, thresholds[thresholds$Group == "CCW", ]) ~ "CCW",
    Species == "Cuculus optatus" & Status == "Breeding" & in_bbox(Longitude, Latitude, thresholds[thresholds$Group == "COW", ]) ~ "COW",
    Species == "Cuculus optatus" & Status == "Breeding" & in_bbox(Longitude, Latitude, thresholds[thresholds$Group == "COE", ]) ~ "COE",
    TRUE ~ "Other"
  ))

# Confirm
breed_assign = ggplot() + 
  geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='grey90',lwd=0.2) +
  coord_sf(xlim = c(-20,185), 
           ylim = c(-30,75), expand = FALSE)+  
  new_scale_fill()+
  new_scale_color()+
  geom_rect(data = thresholds, aes(xmin = longmn, xmax = longmx, ymin=latmn, ymax=latmx, fill=Group),alpha=0.5,col='black',lty=2) +
  geom_point(data=breed_dat %>% filter(Status == 'Breeding'),aes(fill=Group,x=Longitude,y=Latitude,shape=Species),col='black')+
  scale_fill_manual(breaks=c('CCW','CCE','COW','COE','Other'),values=c("#1F78B4","#A6CEE3","#E31A1C","#FB9A99","white" ))+
  scale_shape_manual(values=c(21,24))+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
breed_assign

png('../../figures/ADMIXTURE_StrongholdsK2-Assigned_2024MAY15.png',units='in',res=300,height=4,width=7)
breed_assign
dev.off()

# now that we have assigned each individual a unique ID, let's merge this with the full frame 
breed_id <- breed_dat %>% ungroup %>% select(ID,Group) %>% unique
tracks_df3 <- left_join(tracks_df2,breed_id) %>% 
  #Individuals which did not have june locations should also be assigned 'Other'
  replace_na(list(Group = 'Other'))
tracks_df3 %>% distinct(ID,Group) %>% count(Group)
# # A tibble: 5 × 2
# Group     n
# <chr> <int>
#   1 CCE      20
# 2 CCW     110
# 3 COE       6
# 4 COW       4
# 5 Other   202

#create points 
tracks = st_as_sf(tracks_df3 %>% mutate(Alpha = ifelse(Group == 'Other','Other','Target')), coords = c("Long", "Lat"), 
                  crs = 4326, agr = "constant") 

#create lines 
tracks_lines <- tracks %>% 
  group_by(ID, Species, Group, Alpha) %>% 
  filter(n() >= 2 ) %>%  # Keep only groups with at least two records
  summarize(geometry = st_union(geometry)) %>%
  st_cast("LINESTRING") %>%
  st_as_sf()

migp = ggplot() +
  geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='grey90',lwd=0.2) +
  new_scale_fill()+
  new_scale_color()+
  geom_sf(data=tracks_lines,aes(col=Group,group=ID,alpha=Alpha)) + 
  scale_color_manual(breaks=c('CCW','CCE','COW','COE','Other'),values=c("#1F78B4","#A6CEE3","#E31A1C","#FB9A99","grey40" ))+
  scale_alpha_manual(breaks=c('Other','Target'),values=c(0.2,0.8))+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  coord_sf(xlim = c(-20,185), 
           ylim = c(-30,75), expand = FALSE)+  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
migp

png('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Migration_Map_AssignedK2_2024MAY14.png',units='in',res=300,height=4,width=7)
migp
dev.off()

#### For examining overwintering grounds #####
library(concaveman)
library(smoothr)

#Only select locations in january for our 4 target groups 
winter = tracks_df3 %>% 
  filter(Group != 'Other' & (month(DT) %in% c(1)))
winter_tracks <- st_as_sf(winter, coords = c("Long", "Lat"), crs = 4326, agr = "constant") 

winter_poly = list()
for (group in unique(winter_tracks$Group)) { 
  cat('Making polygon for group: ',group,'\n')
  
  t <- winter_tracks %>% filter(Group == group)
  #concavity is arbitrary, higher number provides a smoother polygon 
  poly <- concaveman(t,concavity=5,length_threshold=2.5) %>% 
    mutate(Group = group)
  winter_poly[[group]] = poly

}

winterp = do.call(rbind,winter_poly)
winter_grounds = ggplot() + 
  geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='grey90',lwd=0.2) +
  new_scale_fill()+
  new_scale_color()+
  geom_sf(data=tracks_lines,aes(col=Group,group=ID,alpha=Alpha)) + 
  scale_color_manual(breaks=c('CCW','CCE','COW','COE','Other'),values=c("#1F78B4","#A6CEE3","#E31A1C","#FB9A99","grey40" ))+
  scale_alpha_manual(breaks=c('Other','Target'),values=c(0.2,0.8))+
  geom_sf(data = winterp,aes(fill=Group),alpha=0.6,col='black',lty=3,lwd=0.5)+
  scale_fill_manual(breaks=c('CCW','CCE','COW','COE','Other'),values=c("#1F78B4","#A6CEE3","#E31A1C","#FB9A99","grey40" ))+
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic(base_size=8)+
  coord_sf(xlim = c(-20,185), 
           ylim = c(-30,75), expand = FALSE)+  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))
winter_grounds

png('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Migration_Map_Overwintering_AssignedK2_2024MAY14.png',units='in',res=300,height=4,width=7)
pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Migration_Map_Overwintering_AssignedK2_2024MAY14.pdf',height=4,width=7)
winter_grounds
dev.off()

# #Calculate distance between mean breeding locations across years
# dist_dat = NULL
# for (id in unique(multi_breed_dat$ID)){
#   s <- multi_breed_dat %>% filter(ID == id)
#   yeardat = NULL
# 
#   #now loop through each consecutive year
#   for (year in seq(1,nrow(s)-1,1)){
#     d <- s %>% mutate(Distance = distHaversine(cbind(Longitude[year],Latitude[year]),cbind(Longitude[year+1],Latitude[year+1]))/1000) %>% ungroup %>% select(Species,ID,Distance) %>% unique
#     yeardat = rbind(yeardat,d %>% mutate(Year = year))
#   }
#   dist_dat = rbind(dist_dat,yeardat)
# }
# 
# #plot the distances observed between sequential years
# dist_dat %>% ggplot(aes(x=Distance,fill=Species))+
#   geom_histogram()+
#   theme_bw()
# 
# #Plot Breeding Locations Over Years for the individuals with low breeding ground fidelity
# plot_breed = tracks_df2 %>%
#   #filter(Status == 'Breeding') %>%
#   filter(Status == 'Overwintering') %>%
#   group_by(ID) %>% mutate(BreedYears = n_distinct(Year_Label)) %>% filter(BreedYears>1)
# plot_ids = dist_dat %>% arrange(desc(Distance)) %>% pull(ID) %>% head(3)
# plot_sub = plot_breed %>% filter(ID %in% plot_ids)
# 
# #Just plot the first 5 individuals
# ggplot() +
#   geom_polygon(data = map.polygon, aes(x = long, y = lat, group = group),fill='white',col='grey90',lwd=0.2) +
#   coord_sf(xlim = c(-20,185),
#            ylim = c(-30,75), expand = FALSE)+
#   new_scale_fill()+
#   new_scale_color()+
#   geom_point(data = plot_sub, aes(x = Long, y = Lat,col=Year_Label,shape=ID), size = 2) +
#   xlab("Longitude") + ylab("Latitude") +
#   scale_shape_manual(values=c(16,3,8))+
#   scale_color_manual(values=viridis(length(unique(plot_sub$Year_Label))))+
#   theme_classic(base_size=8)+
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
#   theme(legend.position = 'top',legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))

```

First, ascribing a bounding box based on the genetic ancestry Q > 99% for west / east.

![ADMIXTURE_StrongholdsK2_2024MAY15](G:\My Drive\Research\Migration\tesselation\2024_05-UnrelatedChyiyin_Fastsimcoal\ADMIXTURE_StrongholdsK2_2024MAY15.png)

Strongholds of > 99% ancestry coefficients for western and eastern groups. Individuals with a mean breeding area within these boxes are assigned to each group. 

![Migration_MonthlyDistances_2024MAY15](G:\My Drive\Research\Migration\tesselation\2024_05-UnrelatedChyiyin_Fastsimcoal\Migration_MonthlyDistances_2024MAY15.png)

Distances traveled in each month

![Migration_MeanMonthlyDistances_2024MAY15](G:\My Drive\Research\Migration\tesselation\2024_05-UnrelatedChyiyin_Fastsimcoal\Migration_MeanMonthlyDistances_2024MAY15.png)

Mean monthly averages between fixes, seems like june and january are the best months to select for breeding and overwintering. 

![ADMIXTURE_StrongholdsK2-Assigned_2024MAY15](G:\My Drive\Research\Migration\tesselation\2024_05-UnrelatedChyiyin_Fastsimcoal\ADMIXTURE_StrongholdsK2-Assigned_2024MAY15.png)

For each individual's tracking data, I calculated the mean lat/long across all years in June. I then intersected this with the bounding boxes to assign each individual into a genetic stronghold. 

We can then color the tracks based on the individual's ancestry:

![Migration_Map_AssignedK2_2024MAY14](G:\My Drive\Research\Migration\tesselation\2024_05-UnrelatedChyiyin_Fastsimcoal\Migration_Map_AssignedK2_2024MAY14.png)

Great, now let's identify overwintering locations for the 4 groups. Unfortunately we do not have ANY january locations for our optatus east samples, and only a single individual for optatus west, so we cannot make polygons for optatus.

Put them together with the migration tracks:

![Migration_Map_Overwintering_AssignedK2_2024MAY14](G:\My Drive\Research\Migration\tesselation\2024_05-UnrelatedChyiyin_Fastsimcoal\Migration_Map_Overwintering_AssignedK2_2024MAY14.png)



## Subset N = 40

```R
### Subset N = 40 based on geography. 
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
all = md %>% filter(Analysis_ADMIXTURE == 1 ) %>% 
  mutate(Species = gsub('C. optatus','Cuculus optatus', 
                        gsub('C. c. canorus','Cuculus canorus',
                             gsub('C. c. bakeri','Cuculus canorus',Species_Latin))),
         LatJit = jitter(Latitude,amount = 1),
         LonJit = jitter(Longitude,amount = 1))

#and assign groups
all = all %>%
  mutate(Group = ifelse(Species == 'Cuculus canorus' & Latitude > 45 & Longitude < 5, 'CanWest',
                        ifelse(Species == 'Cuculus canorus' & Latitude > 35 & Longitude > 135 & Longitude < 165,'CanEast',
                               ifelse(Species == 'Cuculus optatus' & Latitude > 40 & Longitude > 75 & Longitude < 120,'OptWest',
                                      ifelse(Species == 'Cuculus optatus' & Longitude > 135 & Longitude < 150,'OptEast',
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

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Unrelated_Breeding_N190_Distribution_2024APR26.pdf',height=4,width=7)
unrel_breed
dev.off()

all %>% filter(Group != 'Other') %>% group_by(Group) %>% sum_stats(MeanCoverage)

all %>% filter(Group != 'Other') %>% count(Species,Group) %>% 
  ggplot(aes(x=Group,fill=Group,y=n,label=n))+geom_bar(stat='identity')+
  geom_text(vjust=-0.2)+ylab('Count Unrelated Individuals')+xlab('Migratory Group')+
  scale_fill_manual(values=cols)+
  facet_grid(.~Species,scales='free')+
  theme_bw(base_size=18) + theme(legend.position='none')
write.table(all %>% ungroup %>% filter(Group != 'Other') %>%  select(ID,Group),file='~/merondun/cuculus_migration/Samples_Demography_All_CCW-CCE-COW-COE_2024MAY02.pop',quote=F,sep='\t',row.names=F,col.names=F)

#retain n = 10 samples from each stronghold with the highest coverage 
kept = all %>% filter(Group != 'Other') %>% group_by(Species,Group) %>% slice_max(MeanCoverage,n=10)
kept %>% count(Species,Group)
kept %>% group_by(Species,Group) %>% sum_stats(MeanCoverage)
kept %>% count(Group)

write.table(kept %>% ungroup %>% select(ID,Group),file='~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024APR26.pop',quote=F,sep='\t',row.names=F,col.names=F)
write.table(kept %>% ungroup %>% select(ID),file='~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024APR26.list',quote=F,sep='\t',row.names=F,col.names=F)

```

![image-20240426105716120](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20240426105716120.png)



![image-20240426105815823](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20240426105815823.png)

### Subset Samples

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

CHR=$1

mkdir full_vcf

raw_vcf_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/

#From the full gVCF, including invariant sites, subset only samples of interest
bcftools view --threads 10 --types snps --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${raw_vcf_dir}/${CHR}.SNPS.vcf.gz -Ou | \
       bcftools view --min-ac 1 --min-alleles 2 --max-alleles 2 -Oz -o full_vcf/${CHR}.raw.vcf.gz

#Apply filtering, only on MQ, DP
bcftools view --threads 10 --max-alleles 2 -i 'MQ > 30 && INFO/DP > 150 && QUAL > 20' -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz full_vcf/${CHR}.raw.vcf.gz
bcftools index full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz

#Phase VCF with beagle
java -jar -Xmx160g ~/modules/beagle.28Jun21.220.jar gt=full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz out=full_vcf/${CHR}.MQ20-DP150-Q20.PHASED nthreads=8 window=40 overlap=2 impute=true
bcftools index --threads 10 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz

#add the INFO DP,MQ and FMT/DP annotations back onto this VCF, from the pre-phased VCF
bcftools annotate --threads 10 -a full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz -c INFO/DP,INFO/MQ,FMT/DP full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
bcftools index --threads 10 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
```

### ADMIXTURE

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00


#neutral sites
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_breeding/neutral-sites_cuckoo__intergenic-intron-4fold.bed
#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/unrelated_breeding/GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed


MINDP=5
bcftools concat --threads 10 -Ou full_vcf/*Q20.vcf.gz | \
        bcftools sort -Ou | \
        #set genotypes below 5x to missing
        bcftools +setGT -Ou -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #update AC fields after changing genotypes to missing
        bcftools +fill-tags -Ou -- -t AC,AN | \
        bcftools view --min-ac 2 -e 'F_MISSING > 0.1' -Oz -o merged/Autosomes.MQ20-DP150-Q20-MM1-MAC2.vcf.gz

bcftools view --threads 10 merged/Autosomes.MQ20-DP150-Q20-MM1-MAC2.vcf.gz -Ov | \
        #RETAIN neutral sites
        bedtools intersect -header -a - -b $neutral | \
        #EXCLUDE repeats
        bedtools subtract -header -a - -b $repeats | \
        #LD PRUNE
        bcftools +prune -m 0.2 --window 50 -Oz -o merged/Autosomes.MQ20-DP150-Q20-MM1-MAC2-Neutral-LDr2w50.vcf.gz
bcftools index merged/Autosomes.MQ20-DP150-Q20-MM1-MAC2-Neutral-LDr2w50.vcf.gz

~/modules/plink2 --threads 10 --vcf merged/Autosomes.MQ20-DP150-Q20-MM1-MAC2-Neutral-LDr2w50.vcf.gz --chr-set 20 --allow-extra-chr --set-missing-var-ids @:# --make-bed --out merged/Autosomes.MQ20-DP150-Q20-MM1-MAC2-Neutral-LDr2w50

```

### Subset N = 80

Subsample N=20 per group for sensitivity:

```R
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

#and assign groups: OLDER, only N=40
all = all %>%
  mutate(Group = ifelse(Species == 'Cuculus canorus' & Latitude > 45 & Longitude < 5, 'CCW',
                        ifelse(Species == 'Cuculus canorus' & Latitude > 35 & Longitude > 135 & Longitude < 165,'CCE',
                               ifelse(Species == 'Cuculus optatus' & Latitude > 40 & Longitude > 75 & Longitude < 120,'COW',
                                      ifelse(Species == 'Cuculus optatus' & Longitude > 135 & Longitude < 150,'COE',
                                             'Other')))))

#and assign groups: VERY flexible, greater numbers 
all = all %>%
  mutate(Group = ifelse(Species == 'Cuculus canorus' & Latitude > 45 & Longitude < 40, 'CCW',
                        ifelse(Species == 'Cuculus canorus' & Latitude > 35 & Longitude > 100 & Longitude < 165,'CCE',
                               ifelse(Species == 'Cuculus optatus' & Longitude < 120,'COW',
                                      ifelse(Species == 'Cuculus optatus' & Longitude > 120 ,'COE',
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
full_set <- all %>% filter(Group != 'Other') %>% select(ID,Group)
write.table(full_set,file='~/merondun/cuculus_migration/admixture/Population_Designations_Full.pop',quote=F,sep='\t',row.names=F)

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

```

![image-20241023140836565](C:\Users\herit\AppData\Roaming\Typora\typora-user-images\image-20241023140836565.png)

And as before, subset for folded SFS:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=3
#SBATCH --time=200:00:00
# mamba activate snps
# submit as  for i in $(cat Chromosomes.list); do sbatch -J FILTER_${i} ~/merondun/cuculus_migration/snp_calling/branch_folded/7C.Subsample_Demography_N80Subset.sh ${i}; done
CHR=$1

#filtered vcfs directory
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs

outvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n80_vcfs
mkdir $outvcfs $outvcfs/folded $outvcfs/unfolded $outvcfs/invariant_mm1 $outvcfs/raw $outvcfs/stats

#samples
samples=/dss/dsshome1/lxc07/di39dux/merondun/cuculus_migration/deprecated_sampling/double_sampling_2024AUG/Samples_Demography_N20_CCW-CCE-COW-COE_2024AUG12.list

#this bed contains GOOD regions which are neutral for demography
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/extracting_neutral_sites/neutral-sites_cuckoo__intergenic-intron-4fold.bed

#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/GCA_017976375.1_bCucCan1.pri_genomic.CHR_strip.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

### Generate Depth statistics for filtering on DP, retain only sites within mean coverage + 2*sd or mean covearge - 2*sd
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' ${vcfs}/${CHR}_raw.vcf.gz > ${outvcfs}/stats/${CHR}_raw.dp_stats.txt
mean=$(cat ${outvcfs}/stats/${CHR}_raw.dp_stats.txt | datamash mean 3 | xargs printf "%.0f\n")
sd=$(cat ${outvcfs}/stats/${CHR}_raw.dp_stats.txt | datamash sstdev 3 | xargs printf "%.0f\n")
# Round to nearest integer
low=$(echo "$mean - 2*$sd" | bc)
high=$(echo "$mean + 2*$sd" | bc)

#subset samples first, and filter for MQ. Right at the beginning, remove repeats and non-neutral sites
MINDP=5
bcftools view --threads 5 --samples-file ${samples} -Ov ${vcfs}/${CHR}_raw.vcf.gz | \
        #RETAIN neutral sites
        bedtools intersect -header -a - -b $neutral | \
        #EXCLUDE repeats
        bedtools subtract -header -a - -b $repeats | \
        #SET genotypes below MINDP to missing
        bcftools +setGT -Ou -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
        #UPDATE AC fields after changing genotypes to missing
        bcftools +fill-tags -Ou -- -t AC,AN  | \
        #RETAIN only sites with at least 90% genotypes
        bcftools view --threads 5 -i "MQ > 30 & F_MISSING < 0.1 & INFO/DP > ${low} & INFO/DP < ${high}" -Oz -o ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### for FOLDED SFS, just subset
bcftools view --threads 5 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -i 'QUAL > 20' -Oz -o ${outvcfs}/folded/${CHR}_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 ${outvcfs}/folded/${CHR}_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### And then simply grab the invariant sites
bcftools view --threads 5 --max-ac 0 -Oz -o $outvcfs/invariant_mm1/${CHR}_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz
bcftools index --threads 5 $outvcfs/invariant_mm1/${CHR}_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz

### Summarize SNP Counts in tidy format
raw=$(bcftools index -n ${vcfs}/${CHR}_raw.vcf.gz)
neutral=$(bcftools index -n ${outvcfs}/raw/${CHR}_allsitesN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)
folded=$(bcftools index -n ${outvcfs}/folded/${CHR}_snpsN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)
invariant=$(bcftools index -n $outvcfs/invariant_mm1/${CHR}_invariantN80.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz)

echo -e "${CHR}\tRaw\t${raw}" > ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tNeutral_MM1\t${neutral}" >> ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tFolded\t${folded}" >> ${outvcfs}/stats/${CHR}.snp.counts
echo -e "${CHR}\tInvariant\t${invariant}" >> ${outvcfs}/stats/${CHR}.snp.counts

```



## Species Distinction

Subsetting n = 2 from the target and outgroup species, run admixture:

```bash

#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

CHR=chr_1

vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs/unrelated_all/${CHR}_snp.MQ-5X.vcf.gz
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/species_distinction
list=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/species_distinction/ID.list

# Only LD Pruning
bcftools view --threads 10 --samples-file ${list} -Ou ${vcf} | \
        bcftools view --threads 10 --types snps --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -i 'F_MISSING < 0.2' | \
        bcftools +prune -m 0.8 --window 50 -Oz -o ${wd}/${CHR}_snp.MQ-5X-LDr8w50.vcf.gz
bcftools index --threads 10 ${wd}/${CHR}_snp.MQ-5X-LDr8w50.vcf.gz

#Create plink bed files for admixture + perform pca
plink --threads 10 --const-fid --vcf ${wd}/${CHR}_snp.MQ-5X-LDr8w50.vcf.gz --chr-set 29 --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out ${wd}/${CHR}_snp.MQ-5X-LDr8w50

cd ${wd}/admixture

for K in {2..10}; do

echo "RUNNING K: ${K}"

#Run Admixture
admixture -j7 --cv=5 ../${CHR}_snp.MQ-5X-LDr8w50.bed ${K} > ${CHR}_snp.MQ-5X-LDr8w50.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../${CHR}_snp.MQ-5X-LDr8w50 -fname ${CHR}_snp.MQ-5X-LDr8w50.${K}.P -qname ${CHR}_snp.MQ-5X-LDr8w50.${K}.Q -P 10 -o eval_${CHR}_snp.MQ-5X-LDr8w50.${K}

done
```

Also calculate FST:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00


TMPDIR=/dss/dsslegfs01/pr53da/pr53da-dss-0021/tmp
VCF=chr_1_snp.MQ-5X-LDr8w50.vcf.gz

length=1000000000000

tabix ${VCF}
pixy --stats fst --bypass_invariant_check yes --vcf ${VCF} --populations ID.pixypop --window_size ${length} --n_cores 10 --output_folder fst
```

And plot: 

```R
### Plot ADMIXTURE species distinction, n = 2 each 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/species_distinction/all_admix_files')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(stringr)
library(meRo) #devtools::install_github('merondun/meRo')
library(ggpubr)
library(RColorBrewer)

qdir = '.' #directory with Q files
prefixes = gsub('.2.Q','',list.files('.',pattern='.*.2.Q'))
prefix = 'chr_1_snp.MQ-5X-LDr2w50'

counter = 0 

for (prefix in prefixes) { 
  counter = counter + 1
  cat('Working on run: ',prefix,'\n')
    
  admix = melt_admixture(prefix = prefix, qdir = qdir)
  
  #read in metadata
  md = read_tsv('../ID.pixypop',col_names = F)
  names(md) <- c('ID','Species_Latin')
  admixmd = left_join(admix, md %>% select(ID,Species_Latin))
  
  #loop through all admixture runs and extract the average correlation values from evalAdmix, we want to MINIMIZE this! (closest to 0)
  evaldat = NULL; for (Kval in seq(2,10,1)){
    r <- as.matrix(read.table(paste0("eval_",prefix,'.',Kval)))
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
    geom_text(aes(y = 0.15, label = format(signif(median, 2), scientific = TRUE)),size=2) +  
    #ylim(c(-0.015,0.015))+
    geom_hline(yintercept=0,lty=2)+
    geom_point()+ylab('Median +/- IQR Correlation of Residuals') +
    geom_errorbar()+
    theme_bw() + 
    scale_x_continuous(breaks = seq(min(evaldat$K), max(evaldat$K), by = 1)) +
    coord_flip()
  ep
  assign(paste0('e',1),ep)
  
  png(paste0('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20240806_evalAdmix_',prefix,'_SpeciesDistinction.png'),res=300,units='in',height=5,width=6)
  print(e1)
  dev.off()
  
  #Now plot ADMIXTURE 
  admixmd$Species_Latin <- factor(admixmd$Species_Latin,levels=c('C. canorus canorus','C. canorus bakeri','C. optatus','C. saturatus','C. micropterus','C. poliocephalus'))
  adplot =
    admixmd %>% filter(Specified_K <= 9 ) %>%  #specify the levels you want 
    mutate(Specified_K = paste0('K',Specified_K)) %>% 
    ggplot(aes(x = factor(ID), y = Q, fill = factor(K))) +
    geom_col(color = "gray", size = 0.1) +
    facet_grid(Specified_K~Species_Latin, scales = "free", space = "free") +
    theme_minimal(base_size=6) + labs(x = "",y = "") +
    scale_y_continuous(expand = c(0, 0),n.breaks = 3) +
    scale_fill_manual(values=viridis(9,option='turbo'))+
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      #axis.text.x = element_blank(),
      axis.text.x=element_blank(),
      axis.text.y = element_text(size=5),
      panel.grid = element_blank(),
      legend.position = 'bottom',
      plot.title = element_text(size=6)
    )
  adplot
  
  pdf(paste0('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20240806_ADMIXTURE_',prefix,'_SpeciesDistinction.pdf'),height=3,width=3)
  print(adplot)
  dev.off()
  
}
  

### Read in FST
fst <- read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/admixture/species_distinction/FST_Overall.txt')

# Plotting
fst$pop1 <- factor(fst$pop1,levels=c('C. saturatus','C. micropterus','C. poliocephalus','C. canorus bakeri','C. canorus canorus','C. optatus'))
fst$pop2 <- factor(fst$pop2,levels=c('C. saturatus','C. micropterus','C. poliocephalus','C. canorus bakeri','C. canorus canorus','C. optatus'))

fst_plot <- ggplot(fst, aes(x = pop1, y = pop2)) +
  geom_tile(aes(fill = avg_wc_fst), colour = "white") +  # use 'value' because pivot_wider creates it
  geom_text(aes(label = round(avg_wc_fst, 2)),size = 3, vjust = 1) +
  scale_fill_gradient(low = "white", high = "red") +
  facet_grid(~ filter )+
  theme_bw()+
  labs(x = "", y = "", fill = "Avg W&C FST") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fst_plot

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20240806_FST-Heatmap_SpeciesDistinction.pdf',height=3.5,width=8)
fst_plot
dev.off()
```



## MSMC2

### Subset samples

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

CHR=$1

WORKDIR=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop

cd ${WORKDIR}

mkdir full_vcf

raw_vcf_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/

#From the full gVCF, including invariant sites, subset only samples of interest
bcftools view --threads 5 --types snps --samples-file ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list ${raw_vcf_dir}/${CHR}.SNPS.vcf.gz -Ou | \
       bcftools view --min-ac 1 --min-alleles 2 --max-alleles 2 -Oz -o full_vcf/${CHR}.raw.vcf.gz

#Apply filtering, only on MQ, DP
bcftools view --threads 5 --max-alleles 2 -i 'MQ > 30 && INFO/DP > 150 && QUAL > 20' -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz full_vcf/${CHR}.raw.vcf.gz
bcftools index full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz

#Phase VCF with beagle
java -jar -Xmx160g ~/modules/beagle.28Jun21.220.jar gt=full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz out=full_vcf/${CHR}.MQ20-DP150-Q20.PHASED nthreads=8 window=40 overlap=2 impute=true
bcftools index --threads 5 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz

#add the INFO DP,MQ and FMT/DP annotations back onto this VCF, from the pre-phased VCF
bcftools annotate --threads 5 -a full_vcf/${CHR}.MQ20-DP150-Q20.vcf.gz -c INFO/DP,INFO/MQ,FMT/DP full_vcf/${CHR}.MQ20-DP150-Q20.PHASED.vcf.gz -Oz -o full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz
bcftools index --threads 5 full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz

```

Mask for sites below coverage, each a masking file for each individual sample

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

#submit sample as positional for RUN in $(cat ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW-COE_2024MAR13.list); do sbatch -J COV_${RUN} 2.Coverage_Masks.sh ${RUN} ; done
SAMPLE=$1

bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/bams/Illumina_Alignments_Merged
#output mask directory
maskdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/coverage_masks
#directory with the n=40 subsampled phased, neutral and non-repetitive VCFS
vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/full_vcf
#output directory with the individual chromosome level vcfs
indvcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/individual_vcfs

cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop
mkdir $maskdir $maskdir/work $indvcfs

#loop through each chromosome
for CHR in $(cat Chromosomes.list); do

#create a coverage/MQ mask
mosdepth --threads 2 --chrom ${CHR} $maskdir/work/${SAMPLE}_${CHR} $bamdir/${SAMPLE}.bam
avg=$(awk '{print $4}' $maskdir/work/${SAMPLE}_${CHR}.mosdepth.summary.txt | tail -n 1)
min=$(printf "%.0f" $(echo "$avg / 2" | bc -l))
max=$(printf "%.0f" $(echo "$avg * 2" | bc -l))

echo "FOR SAMPLE: ${SAMPLE} average coverage is ${avg}, retainining positions between ${min} and ${max}"
#grab only sites greater than or below half the chromosome avg or double it. grab sites we want to RETAIN!
zcat $maskdir/work/${SAMPLE}_${CHR}.per-base.bed.gz | awk -v x=${max} -v n=${min} '$4 < x && $4 > n' | awk '{OFS="\t"}{print $1, $2, $3}' | bgzip -c > $maskdir/${SAMPLE}_${CHR}.bed.gz

#finally, subset that sample from the population VCF
bcftools view --threads 2 --samples ${SAMPLE} -Ob full_vcf/${CHR}.MQ20-DP150-Q20.PHASED-ANN.vcf.gz | \
        bcftools view --threads 2 --min-alleles 2 --max-alleles 2 --min-ac 1 --max-af 0.9999 -Oz -o $indvcfs/${SAMPLE}_${CHR}.vcf.gz
bcftools index $indvcfs/${SAMPLE}_${CHR}.vcf.gz

done

```

### Run MSMC2

Run MSMC2 and MSMC-IM on all haplotypes (n=8):

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=15
#SBATCH --time=48:00:00

#mamba activate samtools0.1.19
#submit with script 5, otherwise: sbatch ~/merondun/cuculus_migration/msmc/4.Crosscoalescent_Iterative.sh CCW CCE 1

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
#file with the mosdepth coverage masks (sites < 1/2 or > 2 the chromosome-sample-specific coverage are ignored )
maskdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/coverage_masks
#vcfs, individual for sample-chr
vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/individual_vcfs
#genome-wide mappability mask
gwmask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/demography/mappability/masks/
#neutral sites
neutral=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/neutral-sites_cuckoo__intergenic-intron-4fold.bed
#this bed contains BAD sites which are repeats (bed does not include 'low_complexity' and 'simple_repeats', so it's largely TEs)
repeats=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/GCA_017976375.1_bCucCan1.pri_genomic.CHR.AvianRepeatModeler_Repeats-ExcludingLowComplexity-SimpleRepeats.bed

cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop

#submit population 1, population 2, iteration
P1=$1
P2=$2
IT=$3

mkdir -p crosscoal crosscoal/input crosscoal/output

#grab 2 random samples from each population
grep -w ${P1} ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW_20241220.pop | awk '{print $1}' | shuf | head -n 2 > crosscoal/${P1}_${P2}_${IT}.inds
grep -w ${P2} ~/merondun/cuculus_migration/Samples_Demography_N10_CCW-CCE-COW_20241220.pop | awk '{print $1}' | shuf | head -n 2 >> crosscoal/${P1}_${P2}_${IT}.inds

#loop through chromosomes and create the msmc input file
for i in $(cat Chromosomes.list); do

rm crosscoal/${P1}_${P2}_${IT}.${i}_POPmask.list crosscoal/${P1}_${P2}_${IT}.${i}_POPvcf.list

for j in $(cat crosscoal/${P1}_${P2}_${IT}.inds); do

echo "--mask=${maskdir}/${j}_${i}.bed.gz " >> crosscoal/${P1}_${P2}_${IT}.${i}_POPmask.list
echo "${vcfdir}/${j}_${i}.vcf.gz" >> crosscoal/${P1}_${P2}_${IT}.${i}_POPvcf.list

done

generate_multihetsep.py --chr ${i} --negative_mask $repeats --mask ${gwmask}/GCA_017976375.1_${i}.mask.bed.gz --mask $neutral $(cat crosscoal/${P1}_${P2}_${IT}.${i}_POPmask.list) $(cat crosscoal/${P1}_${P2}_${IT}.${i}_POPvcf.list) > crosscoal/input/${P1}_${P2}_${IT}_${i}.multihetsep.txt

done

#Run MSMC2, first on each population, then on both together (exhaustive haplotype comparisons)
~/modules/msmc2_Linux -t 15 -p 1*2+16*1+1*2 -I 0,1,2,3 -o crosscoal/output/${P1}_${P2}_${IT}_msmc_FIRST $(ls crosscoal/input/${P1}_${P2}_${IT}_*)
~/modules/msmc2_Linux -t 15 -p 1*2+16*1+1*2 -I 4,5,6,7 -o crosscoal/output/${P1}_${P2}_${IT}_msmc_SECOND $(ls crosscoal/input/${P1}_${P2}_${IT}_*)
~/modules/msmc2_Linux -t 15 -p 1*2+16*1+1*2 -I  -s -o crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS $(ls crosscoal/input/${P1}_${P2}_${IT}_*)

#Merge all 4 iterations
python ~/modules/msmc-tools/combineCrossCoal.py \
    crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.final.txt \
    crosscoal/output/${P1}_${P2}_${IT}_msmc_FIRST.final.txt \
    crosscoal/output/${P1}_${P2}_${IT}_msmc_SECOND.final.txt > crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.CROSSCOAL_final.txt

#Also run MSMC-IM, sing default settings
mu=1.01e-08
python ~/modules/MSMC-IM/MSMC_IM.py crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.CROSSCOAL_final.txt -p 1*2+16*1+1*2 --printfittingdetails --plotfittingdetails --xlog -mu ${mu} -o crosscoal/output/${P1}_${P2}_${IT}_msmc_ALLHAPS.CROSSCOAL_MSMCIM.final.txt

```

Submit:

```bash
# Loop through each pair of populations
for pair in "CCW CCE" "CCW CO" "CCE CO"; do
  # Split the pair into P1 and P2
  read P1 P2 <<<$(echo $pair)

  # Loop through each iteration
  for IT in {1..10}; do
    # Submit the job with sbatch
    sbatch -J "CROSSCOAL_${P1}_${P2}_${IT}" ~/merondun/cuculus_migration/msmc/4.Crosscoalescent_Iterative.sh ${P1} ${P2} ${IT}
  done
done

```



### Plot

```R
.libPaths('~/mambaforge/envs/R/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/msmc/20241220_CrossCoal3Pop/crosscoal/output')
library(tidyverse)
library(viridis)
library(meRo)
library(RColorBrewer)
library(gmodels)

#for binding later
msmcdat = NULL
files = list.files('.',pattern='ALLHAPS.CROSSCOAL_final.txt')

for (f in files) { 
  cat('Working on file: ',f,'\n')
  t = read_tsv(f)
  string = gsub('_msmc.*','',f)
  splist = strsplit(string,'_')[[1]]
  g1 = split[1]
  g2 = split[2]
  run = split[3]
  t$P1 = g1; t$P2 = g2; t$Iteration = run; t$group = paste0(g1,'__',g2)
  msmcdat = rbind(msmcdat,t)
} 

mu = 1.01e-08
gen = 2.74

#label for y axis 
millions_label <- function(x) {
  return(scales::label_number(suffix = "K", scale = 1e-3)(x))
}

#plot each
msmc_each = msmcdat %>% pivot_longer(!c(Iteration,time_index,left_time_boundary,right_time_boundary,P1, P2, group),names_to = 'pop',values_to = 'lambda') %>% 
  mutate(
    #time = left_time_boundary/mu*gen,
    time = right_time_boundary/mu*gen,
    Ne = (1/lambda)/(2*mu)
  )

msmc_each = msmc_each %>% separate(pop,into=c('d1','ID')) %>% 
  mutate(popID = ifelse(ID == '00',P1,
                        ifelse(ID == '11',P2,'Crosspopulation')),
         popID = gsub('D','G',popID)) %>% 
  dplyr::select(-d1,-ID,-P1,-P2)

#remove time boundaries where the estimate for a single population varies high
unstable = msmc_each %>% group_by(popID,time_index) %>% 
  filter(popID != 'Crosspopulation') %>% 
  summarize(mean=mean(Ne),sd=sd(Ne))
unstp = unstable  %>% 
  ggplot(aes(x=time_index,y=sd,col=popID))+
  geom_point(size=0.75)+
  geom_line()+
  scale_color_viridis(discrete=TRUE)+
  theme_bw(base_size=8)
unstp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_3pop_All-Instability.pdf',height=2,width=3)  
unstp
dev.off()

#inspect Ne for outliers
msmc_each %>% 
  filter(popID != 'Crosspopulation' & time_index >= 3) %>% 
  group_by(group, popID) %>%
  ungroup() %>% 
  ggplot(aes(x=popID,fill=group,y=Ne))+geom_boxplot()+theme_bw()

msmc_each %>% count(group)

# Filter for final plotting, removing the highly unstable time periods (<3)
population_dat = msmc_each %>% filter(popID != 'Crosspopulation' & time_index >= 3)

cols <- data.frame(popID = c('CCW','CCE','CO','COW'),
                   col = c('#33a02c','#b2df8a','gold','#fe9024'))

#plot each 
xp = population_dat %>%
  #filter(popID != 'COW') %>%  
  ggplot(aes(x = time, y = Ne, col = popID,lty=Iteration)) +
  geom_step(lwd = 0.25, show.legend = TRUE) +
  scale_x_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,1e6), labels = c('10KY','100KY','1MY')) + # Modify breaks as needed
  scale_y_log10(limits=c(1e4,3e6),breaks = c(1e4,5e4,1e5,2.5e5,5e5,1e6), labels = millions_label) +  # Modify breaks as needed
  theme(strip.text.y = element_text(angle = 0)) +
  scale_color_manual(values=cols$col,breaks=cols$popID)+
  ylab('Ne')+xlab('Time')+
  theme_test(base_size=7)+theme(legend.position='top')
xp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_3pop_All_NeLines.pdf',height=2,width=2)  
xp
dev.off()

#average Ne of each population
population_dat %>% filter(popID != 'Crosspopulation') %>% group_by(popID) %>% sum_stats(Ne)
population_dat %>% filter(popID != 'Crosspopulation') %>% group_by(popID) %>% slice_min(time)

# popID    mean      sd     se  median    iqr conf_low conf_high
# <chr>   <dbl>   <dbl>  <dbl>   <dbl>  <dbl>    <dbl>     <dbl>
#   1 CCE   105574.  77034.  3596.  90950. 69297.   98508.   112640.
# 2 CCW   190236. 250428. 11689. 124133. 74881.  167265.   213206.
# 3 CO    243492. 363291. 20768. 156906. 93296.  202625.   284359.
# 4 COW   219910. 281023. 16065. 157642. 91971.  188298.   251523.
# # 
# time_index left_time_boundary right_time_boundary Iteration group    lambda   time       Ne popID
# <dbl>              <dbl>               <dbl> <chr>     <chr>     <dbl>  <dbl>    <dbl> <chr>
#   1          3          0.0000292           0.0000494 1         CCW__CCE  187.  13395.  264610. CCE  
# 2          3          0.0000292           0.0000494 1         CCW__CCE   44.0 13395. 1126166. CCW  
# 3          3          0.0000322           0.0000544 1         CCE__CO    29.1 14761. 1698615. CO   
# 4          3          0.0000325           0.0000549 4         CCE__COW   35.1 14898. 1411688. COW  

#or using min/max across runs
popdat_all = msmc_each %>% filter(popID != 'Crosspopulation' & time_index >= 3 & time != Inf) %>% 
  group_by(popID,time_index) %>%
  summarize(time = min(time),
            minNe = min(Ne),
            maxNe = max(Ne),
            meanNe = mean(Ne))

#plot ribbons 
ribplot = popdat_all %>%
  ggplot(aes(x = time, ymin = minNe, ymax=maxNe, fill = popID)) +
  geom_ribbon(alpha=0.5)+
  scale_x_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,1e6), labels = c('10KY','100KY','1MY')) + # Modify breaks as needed
  scale_y_log10(limits=c(1e4,3e6),breaks = c(1e4,5e4,1e5,2.5e5,5e5,1e6), labels = millions_label) +  # Modify breaks as needed
  theme(strip.text.y = element_text(angle = 0)) +
  scale_fill_manual(values=cols$col,breaks=cols$popID)+
  ylab('Min and Max Ne')+xlab('Time')+
  theme_test(base_size=7)
ribplot

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_3pop-All_Ribbons.pdf',height=1.5,width=2.5)  
ribplot
dev.off()


# #interplotated to smooth ribbons 
# n_interpolated_points = 10000
# popdat_all_interpolated <- popdat_all %>%
#   group_by(popID) %>% 
#   mutate(
#     time_new = list(seq(min(time), max(time), length.out = n_interpolated_points)),
#     minNe_new = list(spline(time, minNe, n = n_interpolated_points)$y),
#     maxNe_new = list(spline(time, maxNe, n = n_interpolated_points)$y)
#   ) %>%
#   ungroup() %>%
#   select(-time, -minNe, -maxNe) %>%
#   unnest(c(time_new, minNe_new, maxNe_new)) %>%
#   rename(time = time_new, minNe = minNe_new, maxNe = maxNe_new)
# 
# ribplot = popdat_all_interpolated %>%
#   ggplot(aes(x = time, ymin = minNe, ymax=maxNe, fill = popID)) +
#   #geom_line(data=popdat_all,aes(x=time,y=meanNe))+
#   geom_ribbon(alpha=0.8)+
#   scale_x_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,1e6), labels = c('10KY','100KY','1MY')) + # Modify breaks as needed
#   scale_y_log10(limits=c(1e4,1e6),breaks = c(1e4,1e5,2.5e5,5e5,1e6), labels = millions_label) +  # Modify breaks as needed
#   theme(strip.text.y = element_text(angle = 0)) +
#   scale_fill_viridis(discrete=TRUE)+
#   ylab('Min and Max Ne')+xlab('Time')+
#   theme_test(base_size=7)+
#   theme(legend.position='top')
# ribplot
# 
# pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/figures/MSMC_FourPopulations_Introgression-Error_2024APR1.pdf',height=2.5,width=2.5)  
# ribplot
# dev.off()

#plot cross coalesence rate, and add interspecific vs intraspecific comparison 
msmc = msmcdat %>% 
  #filter(time_index >= 3) %>% 
  #filter(P1 != 'CO' & P2 != 'CO') %>% 
  mutate(time = left_time_boundary/mu*gen,
         crossrate = pmin(1,2*lambda_01 / (lambda_00 + lambda_11))) %>% 
  mutate(Comparison = ifelse(grepl('CC',P1) & grepl('CC',P2),'C. canorus',
                             ifelse(grepl('CO',P1) & grepl('CO',P2),'C. optatus',
                                    'Interspecific')))

#import ice core data 
icecore = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/Ice_Temperature_Reconstructions_Kawamura_NOAA-6076.txt')
icecore = icecore %>% dplyr::rename(time = TopAge)

#scale ice core according to the cross rate
crossrate_range <- range(msmc$crossrate, na.rm = TRUE)
deltaT_range <- range(icecore$deltaT, na.rm = TRUE)
scale_factor <- diff(crossrate_range) / diff(deltaT_range)
icecore$scaled_deltaT <- (icecore$deltaT - mean(deltaT_range)) * scale_factor + mean(crossrate_range)
msmc <- msmc %>% mutate(Opt = ifelse(grepl('COW$',P1) | grepl('COW$',P2),'C. optatus western (n=10)','Merged C. optatus (n=20)'))

climates <- data.frame(
  xmin = c(0, 24e3, 57e3, 122e3),        
  xmax = c(10e3, 34e3, 67e3, 132e3),   
  ymin = -Inf,                           
  ymax = Inf,                            
  fill = c("#c3e6ee","#ffbfbf","#ffbfbf","#c3e6ee"))

#  fill = c("#ffbfbf","#c3e6ee","#c3e6ee", "#ffbfbf"))

ap = msmc %>% 
  #filter(Iteration == 1) %>% 
  ggplot(aes(x = time, y = crossrate,Group=interaction(P1,P2),col=interaction(P1,P2),lty=Iteration)) +
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.1, inherit.aes=FALSE) +
  geom_step(lwd = 0.25, show.legend = TRUE) +
  scale_x_log10(breaks = c(1e4,2.5e4,5e4,1e5,1.5e5,2.5e5), labels = c('10 Ka','25 Ka','50 Ka','100 Ka','150 Ka','250 Ka')) + 
  coord_cartesian(ylim=c(0,1),xlim=c(2e3,3.5e5))+
  scale_y_continuous(
    name = "Crossrate",
    sec.axis = sec_axis(~ . / scale_factor + mean(deltaT_range), name = "Delta T (°C)")
  ) +
  theme(strip.text.y = element_text(angle = 0)) +
  geom_line(data = icecore, aes(x = time, y = scaled_deltaT, group = 1), lwd=0.25,lty=2,inherit.aes=FALSE,color = "black")+
  scale_color_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c'))+
  theme_test(base_size=7)+
  facet_grid(Opt~.)+
  theme(legend.position='top')
ap

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_CrosscoalLines-3pop-All_Migration_IceCore.pdf',height=3,width=3)  
ap
dev.off()

## Rib plot for crosscoal
#or using min/max across runs
popdat_ccr = msmc %>% 
  group_by(time_index,P1,P2,Comparison,Opt) %>%
  mutate(time = min(time)) %>% ungroup %>% 
  group_by(time,P1,P2,Comparison,Opt) %>% 
  sum_stats(crossrate)

#plot ribbons 
ribplot_ccr = popdat_ccr %>%
  ggplot(aes(x = time, ymin = conf_low, Group=interaction(P1,P2), ymax=conf_high, fill = interaction(P1,P2))) +
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.3, inherit.aes=FALSE) +
  geom_ribbon(alpha=0.5)+
  scale_x_log10(breaks = c(1e4,2.5e4,5e4,1e5,1.5e5,2.5e5), labels = c('10 Ka','25 Ka','50 Ka','100 Ka','150 Ka','250 Ka')) + 
  coord_cartesian(ylim=c(0,1),xlim=c(2e3,3.5e5))+
  scale_y_continuous(
    name = "Crossrate",
    sec.axis = sec_axis(~ . / scale_factor + mean(deltaT_range), name = "Delta T (°C)")
  ) +
  scale_fill_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c',"#c3e6ee","#ffbfbf"))+
  theme(strip.text.y = element_text(angle = 0)) +
  geom_line(data = icecore, aes(x = time, y = scaled_deltaT, group = 1),lwd=0.25, lty=2,inherit.aes=FALSE,color = "black")+
  scale_color_manual(values=brewer.pal(3,'Set2'))+
  theme_test(base_size=7)+
  facet_grid(Opt~.)+
  theme(legend.position='top')
ribplot_ccr

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMC_CrosscoalRibbons-3pop-All_Migration_IceCore.pdf',height=2.5,width=3)  
ribplot_ccr
dev.off()


#plot MSMC-IM
imdat = NULL
imfiles = list.files('.',pattern='MSMC_IM.estimates.txt')
#imfiles <- imfiles[!grepl('_CO_',imfiles)]

for (f in imfiles) { 
  cat('Working on file: ',f,'\n')
  t = read_tsv(f)
  string = gsub('_msmc.*','',f)
  split = strsplit(string,'_')[[1]]
  g1 = split[1]
  g2 = split[2]
  run = split[3]
  t$P1 = g1; t$P2 = g2; t$Iteration = run
  imdat = rbind(imdat,t)
} 

mu = 1.01e-08
gen = 2.74
imdat <- imdat %>% mutate(Opt = ifelse(grepl('COW$',P1) | grepl('COW$',P2),'C. optatus western (n=10)','Merged C. optatus (n=20)'))
msmcim = imdat %>% 
  group_by(P1, P2, Iteration,Opt) %>%
  arrange(left_time_boundary) %>%
  mutate(time_index = row_number(),
         time = left_time_boundary*gen) %>%  #note that mu is already integrated into MSMC-IM! 
  #slice(-c(1:4)) %>%
  ungroup() %>% 
  mutate(Comparison = ifelse(grepl('CC',P1) & grepl('CC',P2),'C. canorus',
                             ifelse(grepl('CO',P1) & grepl('CO',P2),'C. optatus',
                                    'Interspecific')))

mp = msmcim %>%  
  ggplot(aes(x=time,y=m,group=interaction(P1,P2,Iteration),lty=Iteration,col=interaction(P1,P2)))+
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.01, inherit.aes=FALSE) +
  geom_step(show.legend = TRUE) +
  scale_color_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c'))+
  #facet_wrap(~Opt,scales='free')+
  scale_x_log10(limits=c(5e3,1e5),breaks = c(5e3,1e4,2e4,5e4,1e5), labels = c('5 Ka','10 Ka','20 Ka','50Ka','100 Ka'))+
  ylab('Migration Rate (m)')+xlab('Time (Years)')+
  theme_test(base_size=7)+
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.2, inherit.aes=FALSE) +
  facet_grid(Opt~.)+
  theme(legend.position='top')
mp

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMCIM_3pop-All_AllLines.pdf',height=2.5,width=3)  
mp
dev.off()


#average by group
mig_dat = msmcim %>% ungroup %>% 
  group_by(P1,P2,Comparison,time_index,Opt) %>%
  mutate(time = min(time)) %>% ungroup %>% 
  group_by(time,P1,P2,Comparison,Opt) %>% 
  sum_stats(m) %>% 
  ungroup %>% 
  arrange(time)

#and just plot 1 for final figure 
mp2 = mig_dat %>%  
  ggplot(aes(x=time,ymin=conf_low,ymax=conf_high,group=interaction(P1,P2),fill=interaction(P1,P2)))+
  geom_ribbon(alpha=0.5)+
  scale_x_log10(limits=c(3000,1e5),breaks = c(5e3,1e4,2e4,5e4,1e5), labels = c('5 Ka','10 Ka','20 Ka','50Ka','100 Ka'))+
  ylab('Migration Rate (m)')+xlab('Time (Years)')+
  theme_test(base_size=7)+
  geom_rect(
    data = climates,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    alpha = 0.3, inherit.aes=FALSE) +
  facet_grid(Opt~.)+
  scale_fill_manual(values=c('darkorchid1','#b2df8a','#33a02c','#b2df8a','#33a02c',"#ffbfbf","#c3e6ee"))+
  theme(legend.position='top')
mp2

pdf('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/figures/20250114_MSMCIM_3pops_All-Ribbons.pdf',height=2.5,width=2.5)  
mp2
dev.off()

```

## LD

```bash

#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00

GROUP=$1

#mamba activate snps
VCF=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/chromosome_vcfs/unrelated_chyiyin/chr_6_snp.MQ-5X-MM1.vcf.gz

plink --vcf $VCF --keep ${GROUP}.plist --double-id --allow-extra-chr \
        --set-missing-var-ids @:# \
        --mac 2 --thin 0.01 --geno 0.1 \
        -r2 gz --ld-window 999999 --ld-window-kb 1000 \
        --ld-window-r2 0 \
        --make-bed --out decay/${GROUP}

zcat decay/${GROUP}.ld.gz | awk '{OFS="\t"}{print $2, $5, $7}' | gzip -c > decay/${GROUP}.ldout.gz

```

```R
### Plot FST Scan 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/fst_scan/all_west_east/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(meRo) #devtools::install_github('merondun/meRo')
library(karyoploteR)

# Load in LD
library(data.table)
ld <- fread('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/fst_scan/all_west_east/LD/decay/W_E.ldout.gz')
lds <- ld %>% filter(BP_A > 28e6 & BP_A < 32e6 & BP_B > 28e6 & BP_B < 32e6)
lds %>% ggplot(aes(x=as.factor(BP_A),y=as.factor(BP_B),fill=R2))+
  geom_tile()+
  theme_void()+
  scale_fill_gradient(low='yellow',high='red')+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  geom_vline(xintercept = c(29.95e6, 31.2e6),col='blue')
```

