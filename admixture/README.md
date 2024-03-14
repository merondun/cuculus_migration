# Population genetic analysis with ADMIXTURE, Tesselation

1. Subset the subgroup (e.g. canorus, or optatus, or can+opt), intersect with neutral (intergenic, 4-fold degenerate), and LD-prune (bcftools +prune -m 0.2 --window 50), create plink bed file.
2. Run ADMIXTURE and evalAdmix. 
3. Plot ADMIXTURE, evalAdmix, and tesselation across the landscape. Can optionally use pie charts. Will require the full `.Q` files for evalAdmix. 

The directory `/full_qs/` contains `.Q` and log files for various subgroup runs. 