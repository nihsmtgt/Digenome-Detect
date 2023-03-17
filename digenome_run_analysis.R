library(tidyverse)
library(doParallel)
cores <- 6
cl <- makeCluster(cores)
registerDoParallel(cl)

setwd("D:/analysis/human_digenome/230317")
output="Sample_DGS_HAP_NT_SIDT1.70x.bed"

NTdata=read_tsv("./digenome-run/Sample_DGS_HAP_NT.70x.bed", col_names = c("chr", "start", "end", "score")) |> 
  mutate(data="RNP(-)")
sampledata=read_tsv("./digenome-run/Sample_DGS_SIDT1.70x.bed", col_names = c("chr", "start", "end", "score")) |> 
  mutate(data="RNP(+)")

data=rbind(NTdata, sampledata)
chr=as_vector(distinct(data, chr))

#delete duplicated detection
#put same ID for duplicated (position difference <10bp) detection
result=foreach(i = chr, .packages = "tidyverse", .combine = "rbind") %dopar% {
  k=1
  temp=data[data$chr==i, ]
  temp=temp[order(temp$start),] |>
    mutate(ID=NULL)
  temp=temp[order(temp$data),]
  temp[1,"ID"] = k
  for (j in 2:nrow(temp)) {
    if(temp[j,"start"]>temp[j-1,"start"]+10){
      k=k+1
    }
    temp[j,"ID"] = k
  }
  temp
}

#keep detection only with highest score
result=result[order(result$score, decreasing = T),]|>
  distinct(chr, data, ID, .keep_all = T)

#delete detection also in the No Treatment data
result=arrange(result, chr, start)
result2=foreach(i = chr, .packages = "tidyverse", .combine = "rbind") %dopar% {
  k=1
  temp=result[result$chr==i, ]
  temp=temp[order(temp$start),] 
  temp[1,"ID"] = k
  for (j in 2:nrow(temp)) {
    if(temp[j,"start"]>temp[j-1,"start"]+10){
      k=k+1
    }
    temp[j,"ID"] = k
  }
  temp
}

result2=dplyr::intersect(arrange(result2, data) |> distinct(chr, ID, .keep_all = T),
                  arrange(result2, desc(data)) |> distinct(chr, ID, .keep_all = T))

write_csv(filter(result2, data=="RNP(+)"), paste0("./", output, "_run_distinct.csv"))
