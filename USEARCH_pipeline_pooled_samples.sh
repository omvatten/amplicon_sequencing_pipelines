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

#----------
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
			
	#The usearch command fastq_mergepairs merges the forward and reverse read for the sample
	#The parameters maxdiffs, minimum merged length, and maximum merged length are specified
	#The merged reads are relabeled with the sample name and saved in a file called "samplename".merged.fastq
	$USEARCH -fastq_mergepairs forw_${samplename}.fq -reverse rev_${samplename}.fq \
	-fastq_maxdiffs ${d} -fastq_pctid ${p} -fastq_minmergelen 200 -fastq_maxmergelen 270 -fastqout ${samplename}.merged.fastq -relabel $samplename.
			
	cat ${samplename}.merged.fastq >> all_d${d}_p${p}.merged.fastq #This adds the merged reads into a file which contains all merged reads from all samples ("all.merged.fastq")
	rm ${samplename}.merged.fastq forw_${samplename}.fq rev_${samplename}.fq #This removed the temporary files
done #The loop goes through all the samples, all merged sequences are placed in the "all.merged.fastq" file and labelled with the sample name

#~ #--- STEP 2: FILTER MERGED READS ---
echo " "
echo "Quality filter"
echo " "
$USEARCH -fastq_filter all_d${d}_p${p}.merged.fastq -fastq_maxee ${e} -fastaout filtered_d${d}_p${p}_e${e}.fa -relabel Filt #This filtering step removes all sequences with an expected error
$USEARCH -fastx_uniques filtered_d${d}_p${p}_e${e}.fa -sizeout -relabel Uniq -fastaout uniques_d${d}_p${p}_e${e}.fa #This picks out the unique sequences from the filtered ones, the abundane of each unique sequence is also noted.

#--- STEP 3: MAKE SVs or OTUs ---
echo " "
echo "Make SVs or OTUs"
echo " "
$USEARCH -unoise3 uniques_d${d}_p${p}_e${e}.fa -zotus zotus_d${d}_p${p}_e${e}_s${s}.fa -minsize ${s}  #This makes SVs
#$USEARCH -cluster_otus uniques_d${d}_p${p}_e${e}.fa -otus otus_d${d}_p${p}_e${e}_s${s}.fa -minsize ${s} -relabel Otu ##Alternatively, this makes OTUs
#$PYT ../scripts/USEARCH_change_Zotu_to_Otu.py #Optionally, this python scripts renames sequences in the "zotus.fa" file so they are called Otu instead of Zotu.

#--- STEP 4: MAKE FREQUENCY TABLE BY MAPPING THE MERGED SEQUENCES FROM THE "ALL.MERGED.FASTQ" FILE TO THE OTUs OR SVs, AND DIVIDING THEM BASED ON SAMPLE WHICH SAMPLE THEY BELONG TO ---
echo " "
echo "Make OTU tables"
echo " "
$USEARCH -otutab all_d${d}_p${p}.merged.fastq -otus zotus_d${d}_p${p}_e${e}_s${s}.fa -id ${id} -otutabout tab_d${d}_p${p}_e${e}_s${s}_i${id}.txt   #Make OTU table based on ZOTUs and all merged sequences

now_end=$(date +"%T") #This records the time the script ended
echo "Run started at: $now_start" #This prints to the terminal the time the scripts started
echo "Run ended at: $now_end" #This prints to the terminal the time the scripts ended
