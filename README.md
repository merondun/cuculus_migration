# Evolutionary Genomics of Migration: Oriental and common cuckoos 

Perpeptual Zenodo DOI: [![DOI](https://zenodo.org/badge/656371438.svg)](https://doi.org/10.5281/zenodo.14873434)

This project aims to investigate the evolutionary genomics underlying migratory strategy in common and oriental cuckoos.

## Project Structure

The repository is organized as follows:

### `/demography/`
Demographic fastsimcoal scripts.

### `/pre_processing/`
Preparing fastq reads, alignment. 

### `/relatedness/`
Identify related individuals. 

### `/snp_calling/`
Calling SNPs from alignments, subsampling individuals. 

### `/species_distinctions/`
Assessing divergence between species using ADMIXTURE and FST.

### `/msmc/`
Run MSMC crosscoalescent and MSMC-IM, using raw gVCFs with basic filtering for invariant sites. 

### `/admixture/`
Run ADMIXTURE and evalAdmix for plotting genetic tesselations across space. 

### `/temperature_periods/`
[Kawamura 2007](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=noaa-icecore-6076) Ice Core Temperature data for identifying warm / cold periods. 

### `/reconstruction/`
Ancestral reconstructions of overwintering and migratory route. 

### `/isolation_by_distance/`
Dividing cuckoos into geographic clusters and estimating FST among these to determine the strength of geographic isolation. 

## Contact

Please open an issue on this repository or reach out Justin Merondun: heritabilities [@] gmail.com, or for demographic modeling Chyi Yin Gwee.
