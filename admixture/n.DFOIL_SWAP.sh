#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00

#mamba activate snps

running_dir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/analyses

#these will be our 4 populations. For sensitivity, also run with P1/2 swapped with P3/4
P1=COW
P2=COE
P3=CCW
P4=CCE
pops=("${P1}" "${P2}" "${P3}" "${P4}" "Outgroup")
names="${P1},${P2},${P3},${P4},Outgroup"
label="${P1}_${P2}_${P3}_${P4}"

#Run 20 iterations for sensitivity 
IT=$1
wd=dfoil/work/${IT}
mkdir dfoil dfoil/out dfoil/work ${wd}

#sample 1 individual from each population
rm $wd/samps_init${IT}.pop
for pop in "${pops[@]}"; do
awk -v p=$pop '$2 == p' SampleSubset_DistanceK_n10_2023DEC06.pop | shuf | head -n 1 >> $wd/samps_init${IT}.pop
awk '{print $1}' $wd/samps_init${IT}.pop > $wd/samps_init${IT}.list
done

#now that we have those individuals we will sample (1 from each pop), subset those from the VCF
#we will rename those individuals simply with the population ID (reheader), create a fasta file 
for CHR in $(cat Chromosomes.list); do
bcftools view --threads 3 --samples-file $wd/samps_init${IT}.list -Ou chromosome_vcfs/${CHR}_snp.MQ-5X-MM1-AA.vcf.gz | \
        bcftools reheader -s $wd/samps_init${IT}.pop | \
        bcftools view --threads 3 --min-alleles 2 --max-alleles 2 --types snps --min-ac 1 --max-af 0.999 -e 'F_MISSING > 0'  -Oz -o $wd/samps_init${IT}.${CHR}.vcf.gz
python ~/modules/vcf2phylip.py --fasta -i $wd/samps_init${IT}.${CHR}.vcf.gz --output-folder $wd
done

#divide into 5000 SNP windows
seqkit concat  $wd/*init${IT}.chr*fasta >  $wd/init${IT}.GW.fa
seqkit sliding -s 5000 -W 5000 -g  $wd/init${IT}.GW.fa -o  $wd/init${IT}.GW5K.fa
mkdir $wd/win_${IT}

# Get unique window ranges from the FASTA headers
windows=$(grep -oP "_sliding:\K\d+-\d+" $wd/init${IT}.GW5K.fa | sort -u)

# Loop through each window range
for window in $windows; do
    # Create a file with headers for the current window
    grep -E ">.*_sliding:$window" $wd/init${IT}.GW5K.fa | sed 's/>//g' > "$wd/win_${IT}/headers_$window.txt"

    # Extract sequences for the current window
    seqkit grep --pattern-file "$wd/win_${IT}/headers_$window.txt" $wd/init${IT}.GW5K.fa > "$wd/win_${IT}/window_$window.fa"
done

rm $wd/win_${IT}/headers*
sed -i 's/_.*//g' $wd/win_${IT}/*fa

#run DFOIL
python3 ~/modules/dfoil/fasta2dfoil.py $wd/win_${IT}/win*fa --out $wd/win_init${IT}.counts --names $names
python3 ~/modules/dfoil/dfoil.py --infile $wd/win_init${IT}.counts --out $wd/win_init${IT}.foil \
    --plot_totals --mincount 2 --plot $wd/win_init${IT}.pdf > $wd/win_init${IT}.precheck
python3 ~/modules/dfoil/dfoil_analyze.py $wd/win_init${IT}.foil > $wd/win_init${IT}.outfoil

# fucntion count warnings in a section
count_warnings() {
    section=$1
    awk "/$section/,/====/"'{print}' "$wd/win_init${IT}.precheck" | grep -c "WARNING"
}

#count warnings for each section
concordant_warnings=$(count_warnings "Checking that concordant patterns")
divergences_warnings=$(count_warnings "Checking that divergences")
terminal_branch_pairs_warnings=$(count_warnings "Checking that terminal branch pairs")

#output
echo -e "${label}\t${IT}\t$concordant_warnings\t$divergences_warnings\t$terminal_branch_pairs_warnings" > dfoil/out/${label}_${IT}.checks

awk -v l=${label} -v i=${IT} '{OFS="\t"}{print $0, l, i}' $wd/win_init${IT}.foil > dfoil/out/${label}_${IT}.foil
