# de-duplication (in-sample and inter-sample) and filtering

NT=read_csv(paste0("./tidy/",output, "_RNP(-).csv"))
RNP=read_csv(paste0("./tidy/",output, "_RNP(+).csv"))

data=rbind(NT, RNP) |>
  filter(CLSCORE>7) |> 
  arrange(desc(CLSCORE)) 

chr = as_vector(distinct(data, chr))

# attach ID when starting positions are close (<10 base)
result = foreach (i = chr,.packages = "tidyverse", .combine = dplyr::bind_rows) %dopar% {
                   k = 1
                   temp = data[data$chr == i,]
                   temp = temp[order(temp$start), ] |>
                     mutate(ID = NULL)
                   temp = temp[order(temp$data), ]
                   temp[1, "ID"] = k
                   for (j in 2:nrow(temp)) {
                     if (temp[j, "start"] > temp[j - 1, "start"] + 10) {
                       k = k + 1
                     }
                     if (temp[j, "data"] != temp[j - 1, "data"]) {
                       k = k + 1
                     }
                     temp[j, "ID"] = k
                   }
                   temp
                 }

# intra-sample duplicated sites
duplicates= result |> 
  group_by(chr,ID) |> 
  filter(n()>1)

# intra-sample de-duplication
result = arrange(result, desc(CLSCORE)) |>
  distinct(chr, data, ID, .keep_all = T)

## inter-sample de-duplication
chr = as_vector(distinct(result, chr))
result = arrange(result, chr, start)
result = foreach(i = chr, .packages = "tidyverse", .combine = dplyr::bind_rows) %dopar% {
                            k = 1
                            temp = result[result$chr == i,]
                            temp = temp[order(temp$start), ]|> 
                              mutate(pair_score=0)
                            temp[1, "ID"] = k
                            for (j in 2:nrow(temp)) {
                              if (temp[j, "start"] > temp[j - 1, "start"] + 10) {
                                k = k + 1
                              }
                              else{
                                temp[j, "pair_score"] = temp[j-1, "CLSCORE"]
                                temp[j-1, "pair_score"] = temp[j, "CLSCORE"]
                              }
                              temp[j, "ID"] = k
                            }
                            temp
                          }

duplicates2= result |> 
  group_by(chr,ID) |> 
  filter(n()>1)

#filtering
result_distinct_filt = mutate(result,
                            filtered = if_else(MQ0 + CLIPS > FwdHead | MQ0 + CLIPS > RevTail, T, F))
RNPf=filter(result_distinct_filt, data=="RNP(+)" & filtered==F) |> 
  arrange(desc(CLSCORE))

write_csv(result_distinct_filt, paste0("./", output, "_distinct_filt.csv"))

