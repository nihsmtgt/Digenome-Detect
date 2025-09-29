library(crisprScore)
library(Biostrings)
library(doParallel)
library(ggplot2)
library(tidyverse)

# preference
output="Sample_DGS_HAP_SIDT1.70x_toolkit2.5" # prefix of output file name
mypath="./input" # path of input files
NTprefix = "_NT" #string which is only in the Negative control file name

sample2_th=2.5 #digenome-toolkit threshold

#reference genome
refGenome=readDNAStringSet("~/workspace/hg38/hg38.fa.bgz") # read reference genome

#genome editing tool
PAM = DNAString("NGG")
sgRNAseq = DNAString("GTCTGCTGGCGAACCACAACA")

ToolkitComparison = T # If True, Digneome-detect data will be compared to Digenome-toolkit data
if(ToolkitComparison == T){
  NTdataname="./digenome-run/Sample_DGS_HAP_NT_digenome-seq_result.0.01.bed"
  sampledataname="./digenome-run/Sample_DGS_SIDT1.70x.0.bam_digenome-seq_result.0.01.bed" #file path of the each data of Digenome-toolkit
}

doGGGenome = T # If True, use GGGenome data for distance
if(doGGGenome ==T){
  gg = read_csv("./SIDT1_d6_new.csv.gz") #GGGenome search result
}

cores <- 24 # No. of thread used

#main
cl <- makeCluster(cores)
registerDoParallel(cl)

#digenome-toolkit
if(ToolkitComparison == T){
  source("./1_digenome_run_analysis_01.R", echo = T)
}

#digenome-detect
source("./1_BedtoCsv.R", echo = T)
source("./2_deduplication.R", echo = T)
output=paste0(output, "_distinct_filt")
source("./3_threshold.R", echo = T)

#comparison
if(ToolkitComparison == T){
  source("./4_compare.R", echo = T)
  output=paste0(output, "_compared")
}

#get sequence and distance
source("./5_getsequences.R", echo = T)
output=paste0(output, "_seq")
if(doGGGenome == T){
  source("./6_referGGGenome.R", echo = T)
  output=paste0(output, "_gg")
}
source("7_alignment2.R", echo = T)

stopCluster(cl)

