genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa

#create dict
gatk CreateSequenceDictionary -R $genome
#split genome by Ns 
gatk ScatterIntervalsByNs -R $genome -O scatter.interval_list -OT ACGT

#split intervals
gatk SplitIntervals -R $genome -L scatter.interval_list --scatter-count 100 -O interval_files/
