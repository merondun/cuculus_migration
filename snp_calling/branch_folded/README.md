# Re-filter for the folded SFS

7B. `7B.Subsample_Demography_N40Subset.sh`, re-filters raw VCFs with a slightly different logic, as follows:

- a. retain only neutral sites - either intergenic or 4-fold degenerate
- b. remove repeats (largely TEs, keep simple repeats)
- c. set genotypes below 5x to missing for each sample
- d. retain sites with `MQ > 30` and `F_MISSING < 0.1`, ensuring at least 90% scored genotypes. Also remove sites based on DP: Calculate mean and sd coverage for each chromosome, and remove sites `> mean + 2*sd` or `< mean - 2*sd` 
- e. now subset folded SNPs
- f. for unfolded, assign derived allele based on genotypes in the outgroups. Only sites where all outgroups have the same allele are assigned as derived. Set the derived allele to the REF field, and remove any unassigned alleles
- g. also grab invariant sites. This file will have less SNPs than in the directory above, because the <5X genotype was also applied here to invariant sites - whereas it wasn't in the previous run (only the 10% missingness filter). 

8B. `8B.Merge_Autosomes.sh`, merges chromosome-level VCFs into whole-autosomal files for folded, unfolded, and invariant sites.

9B. `9B.SubsetOnlyN40.sh`, this script just subsets only the N=40 samples, because the outgroup species are present in the previous files. This will add the sites which are invariant within the N=40 into the invariant sites file, and create folded and unfolded autosomal files. **It will also create an additional MAF5% and prune for LD for admixture analyses, for both folded and unfolded**. 


## Full file paths

**SNP file for folded SFS: N = 29,076,495 SNPs**

`/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz`

**SNP file for unfolded SFS: N = 14,032,074 SNPs**

`/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA.vcf.gz`

**Invariant sites: N = 691,881,386 invariant sites**

`/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs/merged/n40/Autosomes_invariantN40.MQ-5X-MM1-Neutral-NonRepetitive.vcf.gz`

## LD-Pruned files for ADMIXTURE 

**SNP file, equivalent to folded SFS filter except with LD pruning: N = 2,004,430 SNPs** 

`/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs/merged/n40/Autosomes_Folded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-MAF5-LDr2w50.vcf.gz` 

**SNP file, equivalent to unfolded SFS filter with LD pruning: N = 982,381 SNPs**

`/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses/subsampled_n40_vcfs/merged/n40/Autosomes_Unfolded_snpsN40.MQ-5X-MM1-Neutral-NonRepetitive-AA-MAF5-LDr2w50.vcf.gz` 