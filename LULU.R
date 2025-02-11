# LULU script

#uses LULU to clean OTU tables 
library(dplyr)
library(lulu)
library(compiler)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3){stop("Not enough commandline args!")}

matchL = args[1]
otuF = args[2]
fastaF = args[3]
logD = args[4]

matches = read.table(matchL,header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
otuM = read.table(otuF, header=TRUE, comment.char = '&', sep="\t", row.names=1, as.is=TRUE)

lulu_results <- lulu(otuM, matches, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)

write.table(lulu_results$curated_table,quote =FALSE,sep="\t",file=paste0(logD, "/lulu_curated_table.txt", sep=""),col.names = NA)
write.table(lulu_results$curated_otus, quote =FALSE, sep="\t", file=paste0(logD, "/lulu_curated_otus.txt", sep=""), col.names = NA)
write.table(lulu_results$discarded_otus,quote =FALSE,sep="\t",file=paste0(logD, "/lulu_discarded_otus.txt", sep=""),col.names =FALSE,row.names =FALSE)

print("Done.")
