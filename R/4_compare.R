
sample1=read_csv(paste0(output, ".csv"))
sample1name="detect"
sample1_th=read_csv(paste0(output, "_threshold.csv"))$threshold_low
sample2=read_csv(toolkitfilename)
sample2name="toolkit"

sample1=sample1|> 
  filter(CLSCORE>sample1_th & filtered==F & pair_score<=7) |> 
  mutate(score=CLSCORE) |> 
  select(c(any_of(names(sample2)), "score")) 

sample2=sample2|> 
  select(any_of(names(sample1))) |> 
  filter(score>sample2_th & pair_score<=sample2_th)

compare= rbind(mutate(sample1, group="sample1"),
               mutate(sample2, group="sample2"))

#start, endをそろえる
chr=as_vector(distinct(compare, chr))
k=0
for (i in chr) {
  k=k+1
  temp=compare[compare$chr==i, ]
  temp=temp[order(temp$start),]
  temp[1, "ID"]=k
  for (j in 2:length(temp$start)) {
    if(temp[j,"start"]>temp[j-1,"start"]+10){
      k=k+1
    }
    temp[j, "ID"]=k
  }
  if (i == chr[1]){
    result=temp
  }
  else{
    result=rbind(result, temp)
  }
}

for (i in 1:max(result$ID)) {
  if(nrow(filter(result, ID==i))!=1){
    result[result$ID==i, "start"]= min(result[result$ID==i, "start"])
    result[result$ID==i, "end"]= max(result[result$ID==i, "end"])  
  }
}

#まとめる
compare2=pivot_wider(result, names_from = "group", values_from = "score") |> 
  replace_na(list(sample1=0, sample2=0)) 

compare3=mutate(compare2, TFsample1=sample1>sample1_th, TFsample2=sample2>sample2_th) |> 
  dplyr::rename(!!sample1name:=sample1, !!sample2name:=sample2, !!paste0("TF",sample1name):=TFsample1,!!paste0("TF",sample2name):=TFsample2)

write_csv(compare3, paste0("./",output,"_compared.csv"))
