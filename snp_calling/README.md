# Basic scripts for calling SNPs

1. Divide genome into blocks 
2. Call SNPs on each block
3. Merge each block by chromosome. 
4. Filter chromosomes, including all samples.
5. Summarize SNPs
6. Merge VCFs into an all sample autosomal file.
7. Subsample the n=40 samples for demographic inference.
8. Merge subsampled VCFs into a single autosomal file. 

The file `Outgroups.list` is a list of samples used for polarization. 