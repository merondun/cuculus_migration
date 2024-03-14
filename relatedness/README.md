# Removing related individuals

1. Identify relatives from plink using MAF 5% SNPs on chromosome 10 with a threshold of 0.185, following 10.1038/nprot.2010.116
2. Remove individuals based on missingness: if sample A is missing more data (calculated from vcftools), then remove A, otherwise B.

The remaining list files are simply lists of samples by each subgrouping. `.plink` files can be used for subsetting files with plink. 