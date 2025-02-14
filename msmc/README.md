# Coalescent analyses

`1.Mappability_Mask.sh` Prep genome, create mappability mask

`1B.Mappability_Mask.py` Python script required for mappability mask

`2.Subset_Samples_Phase_VCF.sh` Using the all-sites raw VCF, subset samples, filter, and phase

`3.Coverage_Mask_Sample_Level.sh` Create a sample and chromosome-specific coverage mask file (with retained sites)

`4.Crosscoalescent_Iterative.sh` Run crosscoalescent and MSMC-IM, iteratively

`5.Submit_MSMC2.sh` Submit script 4 using all the population pairs across iterations

`6.PlotResults.R` Plot the results, including MSMC2 and MSMC-IM

All required output files are in `/outputs/`

