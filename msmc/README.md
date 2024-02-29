# Coalescent analyses

1. Prep genome, create mappability mask
1A. Python script required for mappability mask
2. Using the all-sites raw VCF, subset samples, filter, and phase
3. Create a sample and chromosome-specific coverage mask file (with retained sites)
4. Run crosscoalescent and MSMC-IM, iteratively