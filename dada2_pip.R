# Script to run DADA2
# based on the DADA2 tutorial (https://benjjneb.github.io/dada2/tutorial.html)

#packageVersion("dada2")
#library(reticulate)
#conda_list()

library(dada2)
library(dplyr)
library(stringr)
library(seqinr)

#ARGS parsing
args = commandArgs(trailingOnly=TRUE)

path_output = args[1]
omega = as.numeric(args[2])
ncores = as.integer(args[3])
bp4error = 5e7

print("DADA2 is executed with OMEGA_A: ")
print(omega)

# has to be the dataset before dereplication
pathF <- "/Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/LENGTH416_420"

# get all file names
path_files <- list.files(pathF, pattern=".fq", full.names=TRUE)
sample_names <- sapply(strsplit(basename(path_files), "_"), `[`, 1)
names(path_files) <- sample_names

# remove Ns from fastqs
filtered_path_R1 <- file.path("/Datos/mockcommunities_EPA45/epa45_metabarcodes/sequences/02_preprocessed/preprocess_fasta_for_dada2/", paste0(sample_names, "_filtered.fastq.gz"))
names(filtered_path_R1) <- sample_names

# dereplicate the reads
if (!dir.exists(path_output)) {dir.create(path_output)}
derep_run1 <- derepFastq(filtered_path_R1, verbose=TRUE)
names(derep_run1) <- sample_names

# calculate error model (dada2 takes error models into account)
try_error <- try(err1 <- learnErrors(derep_run1, nbases = bp4error, multithread=ncores))
if (class(try_error) == "try-error") {err1 <- learnErrors(derep_run1, nbases = bp4error, multithread=1)}

# run dada2
# default for omega_a is 1e-40
dadaResult_run1 <- dada(derep_run1, err=err1, multithread=ncores, pool=FALSE, OMEGA_A=omega)

# create Sequence table
seqtab <- makeSequenceTable(dadaResult_run1)

# write fasta
seqnum <- paste0("ASV_", seq(ncol(seqtab)))
uniqueSeqs <- as.list(colnames(seqtab))
tmp_name <- paste("dada2_", as.character(omega), ".fasta", sep="")
fasta_file <- paste(path_output, tmp_name, sep="")
write.fasta(uniqueSeqs, seqnum, fasta_file)
seqtab_path <- paste(path_output, "seqtab.csv", sep="")
write.table(seqtab, seqtab_path, sep=",")

# Flip table
seqtab.t_default <- as.data.frame(t(seqtab))
rep_set_ASVs <- as.data.frame(rownames(seqtab.t_default))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t_default)` 
rep_set_ASVs$`rownames(seqtab.t_default)` <- NULL


# Add ASV numbers to table
rownames(seqtab.t_default) <- rep_set_ASVs$ASV_ID
head(seqtab.t_default)

# write abundance table
tmp_name <- paste("abundance_table_", as.character(omega), ".csv", sep="")
abundance_table_path_default <- paste(path_output, tmp_name, sep="")
write.csv(seqtab.t_default, abundance_table_path_default, sep=",", quote=FALSE)

print("done.")
