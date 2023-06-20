#!/bin/bash

#In a directory with the cleaned .fq.gz files, this will create ID (for us just SRR ID), SM (most important, for VCF files), LB (in case same library across multiple lanes), and PU (run data) fields. It will also count the number of SRRs per SM for merging later. 

mkdir RGs

for i in $( ls *.trim.fastq.gz | sed 's/\..*//g' ); do
        zcat $i.trim.fastq.gz | head -n 1 > ./RGs/$i.txt
        awk 'BEGIN {FS=":"}{print FILENAME, $3, $4, $NF}' ./RGs/$i.txt > ./RGs/$i.tmp
        ex -sc '%s/.txt//g' -c x ./RGs/$i.tmp
        ex -sc '%s:./RGs/::g' -c x ./RGs/$i.tmp

        #okay, now subet the parts we want
        sed 's/__.*//g' ./RGs/$i.tmp > ./RGs/$i.SM
        sed 's/.*__//g' ./RGs/$i.tmp | cut -d ' ' -f 1 > ./RGs/$i.ID
        awk '{print $4}' ./RGs/$i.tmp > ./RGs/$i.tmp2
        paste ./RGs/$i.SM ./RGs/$i.tmp2 | sed 's/\t/./g' > ./RGs/$i.LB

        #run details for PU
        zcat $i.trim.fastq.gz |head -n 1 |cut -f 3-5 -d ":" > ./RGs/$i.PU

        #Determine number of runs per sample, and output to .numlibs file
        ID="$(cat ./RGs/${i}.ID | tail -1)";
        SM="$(cat ./RGs/${i}.SM | tail -1)";
        echo ${ID} >> ./RGs/${SM}.libco
        cat ./RGs/${SM}.libco | sort -u > ./RGs/${SM}.numlibs

done

rm ./RGs/*tmp*
rm ./RGs/*txt*
