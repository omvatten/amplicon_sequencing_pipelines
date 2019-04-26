#R script implementing DADA2
library(dada2)

start_date <- date()
print(start_date)

#Path to fastq files
path <- "~/path_to_folder_with_sequence_data"
list.files(path)

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#Look at the quality score for forward and reverse read, two samples are taken as example for each
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in dada_filtered/ subdirectory
filt_path <- file.path(path, "dada_filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#The truncLen values are set based on the quality plots above (the position for forw and rev where quality drops)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,150),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out) #Check how many reads passed the filter

#The algorithm needs to learn the error rates in these specific samples
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Run dada algorithm
# We tested varying OMEGA_A (1e-10 or 1e-40) and pool (TRUE or FALSE)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, OMEGA_A=1e-40, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, OMEGA_A=1e-40, pool=TRUE)
dadaFs[[1]]

#Merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Write frequency table to text file
write.table(seqtab.nochim,"~/path/DADA2_frequency_table.txt",sep="\t")

# Print start time and end time
end_date <- date()
print(start_date)
print(end_date)
