# Population genetic analysis with ADMIXTURE & Tesselation on N=40 Samples

The script `2.ADMIXTURE_N40.sh` just runs ADMIXTURE and evalAdmix on the subsetted N=40 samples. 

It runs K2-10 on 2 subsets: a larger SNP set which was used also for the folded spectra, and a smaller SNP set which only includes alleles that have a known derived state, based on homozygous state in the outgroup C. poliocephalus and C. micropterus. 

The directory `/full_qs/` contains `.Q` and log files for the 2 runs. 