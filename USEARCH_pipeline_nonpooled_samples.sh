#!/bin/bash

#The script is run in bash, the line above is needed for bash scripts.
#To run the script, in the terminal go to the fold where the script is located, type bash "script name", then enter

USEARCH=usearch #path to usearch
PYT=python3  #path to python3
DATA=../../fastq_data #path to directory where fastq data files are located

rm -rf ../out #If there is a fold called "out" in the directory, this is removed
mkdir -p ../out #This creates a new, empty folder called "out"
cd ../out #This means were are operating in the out folder (i.e. the generated files will be placed there, unless otherwise specified)

now_start=$(date +"%T") #Records the time the script started
echo "Current time start: $now_start" #Prints the time to the terminal

d=50 #fastq_maxdiffs parameter value
p=80 #fastq_pctid parameter value
e=1 #fastq_maxee parameter value
s=8 #minsize parameter value
id=0.99 #id parameter value

echo " "
echo "Create sample list"
echo " "
readarray sample_names_list < ../scripts/sample_sheet #Reads names of samples from a file called "sample_sheet" (located in the scripts folder) into an array called sample_names
endforw="_R1_001.fastq" #Specifies the ending of the forward fastq file names
endrev="_R2_001.fastq" #Specifies the ending of the reverse fastq file names

#--- STEP 1: MERGE READS ---
echo " "
echo "Merge reads"
echo " "
	
for samplename in ${sample_names_list[@]} #This loop goes through the samples one by one
do
	cat $DATA/${samplename}$endforw >forw_${samplename}.fq #Save the fastq file with the forward reads in a temporary file names forw_"samplename".fq
	cat $DATA/${samplename}$endrev >rev_${samplename}.fq #Save the fastq file with the reverse reads in a temporary file names rev_"samplename".fq
	
	$USEARCH -fastq_mergepairs forw_${samplename}.fq -reverse rev_${samplename}.fq \
	-fastq_maxdiffs ${d} -fastq_pctid 70 -fastq_minmergelen 200 -fastq_maxmergelen 270 -fastqout merged.${samplename}.fastq -relabel $samplename.


	#~ #--- STEP 2: FILTER MERGED READS ---
	echo " "
	echo "Quality filter"
	echo " "
	$USEARCH -fastq_filter merged.${samplename}.fastq -fastq_maxee ${e} -fastaout filtered.${samplename}_d${d}_e${e}.fa -relabel Filt #This filtering step removes all sequences with an expected error exceeding e
	$USEARCH -fastx_uniques filtered.${samplename}_d${d}_e${e}.fa -sizeout -relabel Uniq -fastaout uniques.${samplename}_d${d}_e${e}.fa #This picks out the unique sequences from the filtered ones, the abundane of each unique sequence is also noted.

	#--- STEP 3: MAKE ZOTUS USING THE UNOISE ALGORITHM ---
	echo " "
	echo "Make SVs or OTUs"
	echo " "
	$USEARCH -unoise3 uniques.${samplename}_d${d}_e${e}.fa -zotus zotus.${samplename}_d${d}_e${e}_s${s}.fa -minsize ${s}
	#$USEARCH -cluster_otus uniques.${samplename}_d${d}_e${e}.fa -otus otus${samplename}_d${d}_e${e}_s${s}.fa -minsize ${s} -relabel Otu ##Alternatively, this makes OTUs
	#$PYT ../scripts/USEARCH_change_Zotu_to_Otu.py #Optionally, this python scripts renames sequences in the "zotus.fa" file so they are called Otu instead of Zotu.

	#--- STEP 4: MAKE AND ZOTU TABLE BY MAPPING THE MERGED SEQUENCES FROM THE "ALL.MERGED.FASTQ" FILE TO THE ZOTUS, AND DIVIDING THEM BASED ON SAMPLE WHICH SAMPLE THEY BELONG TO ---
	echo " "
	echo "Make OTU tables"
	echo " "
	$USEARCH -otutab merged.${samplename}.fastq -otus zotus.${samplename}_d${d}_e${e}_s${s}.fa -id ${id} -otutabout tab.${samplename}_d${d}_e${e}_s${s}_i${id}.txt   #Make OTU table based on ZOTUs and all merged sequences
done			

$PYT ../scripts/USEARCH_merge_nonpooled_tables.py

now_end=$(date +"%T") #This records the time the script ended
echo "Run started at: $now_start" #This prints to the terminal the time the scripts started
echo "Run ended at: $now_end" #This prints to the terminal the time the scripts ended

