data=read_csv(paste0(output, ".csv"))

data2=mutate(data, seq="N") 

for (i in 1:nrow(data)) {
  temp=data[i, ]
  data2[i, "seq"]=as.character(subseq(refGenome[temp$chr], temp$start-30, temp$end+30))
}

write_csv(data2, paste0("./", output, "_seq.csv"))
