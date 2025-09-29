
#設定
data=read_csv(paste0(output, ".csv"))

gg2=arrange(gg, distance, desc(PAMmatch)) |> 
  select(c("name",
           "position",
           "position_end",
           "distance",
           "ins",
           "del",
           "PAMmatch"))|>
  dplyr::rename(ggstart = position, ggend = position_end)

notmatch=tibble(
  ggstart = -1,
  ggend = -1,
  distance = -1,
  ins = -1,
  del = -1,
  PAMmatch = FALSE,
  PAM_dist = -1
)

result = foreach (i = 1:nrow(data),.packages = "tidyverse",.combine = dplyr::bind_rows) %dopar% {
                           temp2 = data[i, ]
                           temp = filter(gg2,
                                         ggstart <= temp2$start & ggend >= temp2$end & name == temp2$chr)
                           if (nrow(temp) != 0) {
                             temp_r = temp[1, -1] |>
                               mutate(PAM_dist = min(temp[temp$PAMmatch == TRUE, ]$distance))
                           } else{
                             temp_r = notmatch
                           }
                           temp_r
                         }

result_bind = cbind(data, result)
result_bind[result_bind$PAM_dist==Inf,"PAM_dist"]=-1
write_csv(result_bind, paste0("./", output, "_gg.csv"))

