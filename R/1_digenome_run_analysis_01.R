
NTdata = read_tsv(
  NTdataname,
  col_names = c("chr", "start", "end", "score")
) |>
  mutate(data = "RNP(-)")  |> 
  mutate(chr = str_replace(chr, "chrchr", "chr")) #this data contained excess strings in chr column

sampledata = read_tsv(
  sampledataname,
  col_names = c("chr", "start", "end", "score")
) |>
  mutate(data = "RNP(+)")

data = rbind(NTdata, sampledata) |> 
  filter(score>sample2_th)
chr = as_vector(distinct(data, chr))

#delete duplicated detection
#name same ID for duplicated (position difference <10bp) detection
result = foreach(i = chr,
                 .packages = "tidyverse",
                 .combine = "rbind") %dopar% {
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

#keep detection only with highest score
result = result[order(result$score, decreasing = T), ] |>
  distinct(chr, data, ID, .keep_all = T)

#delete detection also in the No Treatment data
result = arrange(result, chr, start)
result2 = foreach(i = chr,
                  .packages = "tidyverse",
                  .combine = "rbind") %dopar% {
                    k = 1
                    temp = result[result$chr == i,]
                    temp = temp[order(temp$start), ] |> 
                      mutate(pair_score=0)
                    temp[1, "ID"] = k
                    
                    for (j in 2:nrow(temp)) {
                      if (temp[j, "start"] > temp[j - 1, "start"] + 10) {
                        k = k + 1
                      }
                      else{
                        temp[j, "pair_score"] = temp[j-1, "score"]
                        temp[j-1, "pair_score"] = temp[j, "score"]
                      }
                      temp[j, "ID"] = k
                    }
                    temp
                  }

write_csv(result2,
          paste0("./", output, "_toolkit_distinct.csv"))
toolkitfilename=paste0("./", output, "_toolkit_distinct.csv")

