# Evolutionary Genomics of Migration: Oriental and common cuckoos 

This project aims to investigate the evolutionary genomics underlying migratory strategy in common and oriental cuckoos.

## Project Structure

The repository is organized as follows:

### `/pre_processing/`: Preparing fastq reads, alignment. 

### `/snp_calling/`: Calling SNPs from alignments, subsampling individuals. 

### `/relatedness/`: Identify related individuals. 

### `/msmc/`: Run MSMC crosscoalescent and MSMC-IM, using raw gVCFs with basic filtering for invariant sites. 

### `/admixture/`: Run ADMIXTURE and evalAdmix for plotting genetic tesselations across space. 

### `/pi_tajima/`: Estimating π and Tajima's D on the n = 40 subsetted samples used for demographic inference.  

### `/temperature_periods`: [Kawamura 2007](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=noaa-icecore-6076) Ice Core Temperature data for identifying warm / cold periods. 

## Contact

For any questions or concerns, please open an issue on this repository or reach out to me (Justin Merondun): heritabilities@gmail.com
