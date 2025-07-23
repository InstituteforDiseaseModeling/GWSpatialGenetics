#!/bin/bash

#Worm pipeline name manifest generator
#CORNELL only

#mimic the thing Jessica wrote

#desired output columns
#fq1, fq2, sample,unit,ena_accession

#fq1 - forward read location
#fq2 - reverse read location
#sample - example: CHDDOG2016_01148 
## in this case, between 4th and 6th _
#unit - just put 1
#ena_accession - copy sample, add '-1' after it
date
echo -e "Creating metadata.tsv for Cornell style mtDNA GW samples.\nWriting header."
#write header for metadata.tsv
echo -e "fq1\tfq2\tsample\tunit\tena_accession" > metadata.tsv

#get working directory
workdir=$(pwd|sed 's/$/\//')
echo "Adding individual samples"
#iterate through list of files
filenames=$(ls *.fastq.gz |cut -d"_" -f1-8| sort | uniq)
while IFS= read -r a; do ##read list of file IDs one at a time
	#output all columns per unique sample [e.g. not R1 or R2, but both]
	#get 'sample'
	samp=$(echo $a|cut -d"_" -f5-6)
	echo -e $workdir$a"_R1.fastq.gz\t"$workdir$a"_R2.fastq.gz\t"$samp"\t1\t"$samp"-1" >> metadata.tsv
	echo "Processed sample "$samp"."
done <<< "$filenames" 

echo "Done at `date`"
