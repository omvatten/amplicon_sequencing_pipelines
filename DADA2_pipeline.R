#R script implementing DADA2
library(dada2)

start_date <- date()
print(start_date)

#Ensure files are unzipped ("gzip -d *.gz" in commandline)
path <- "~/path_to_folder_with_sequence_data"
list.files(path)

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))
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

mergers <- vector("list", length(sample.names))
dadaF <- vector("list", length(sample.names))
dadaR <- vector("list", length(sample.names))
names(mergers) <- sample.names
names(dadaF) <- sample.names
names(dadaR) <- sample.names

for(sam in sample.names){
  cat("Processing:", sam, "\n")
  
  #Dereplication
  derepF <- derepFastq(filtFs[[sam]], verbose=TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose=TRUE)
  dadaF[[sam]] <- dada(derepF, err=errF, multithread=TRUE, OMEGA_A=1e-40)
  dadaR[[sam]] <- dada(derepR, err=errR, multithread=TRUE, OMEGA_A=1e-40)
  mergers[[sam]] <- mergePairs(dadaF[[sam]], derepF, dadaR[[sam]], derepR, maxMismatch=0, verbose=TRUE)
}

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))

#Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track

#Save as text file
write.table(seqtab.nochim,"~/path/DADA2_frequency_table.txt",sep="\t")

end_date <- date()
print(start_date)
print(end_date)

