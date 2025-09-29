# convert separated bed files to single csv for R programs

# main
files=list.files(path=mypath, full.names = TRUE)

result=foreach(i = 1:length(files), .packages = "tidyverse", .combine = dplyr::bind_rows) %dopar% {
  input=files[i]
  bed=read_tsv(input, col_names = c("chr", "start", "end", "INFO"))
  
  sep=letters 
  sep=sep[1:str_count(bed[1, "INFO"], "=")]
  
  bed2=bed |> 
    mutate(INFO=str_replace_all(INFO, ":", ";"))|>
    separate(INFO, sep =";" , into=sep) 
  bed2=bed2 |>
    pivot_longer(cols=all_of(sep), names_to = "name", values_to = "value")|>
    separate(value, into=c("name", "value"), sep = "=", convert=T)|>
    pivot_wider(names_from = "name", values_from = "value")
  filename=last(last(str_split(input, "/")))
  bed2=mutate(bed2, data=if_else(str_detect(filename, NTprefix), "RNP(-)", "RNP(+)"))
}

if(!dir.exists("./tidy")){
  dir.create("./tidy")
}


for (i in 1:nrow(distinct(result, data))) {
  write_csv(result[result$data==as.character(distinct(result, data)[i,]),], paste0("./tidy/",output, "_", distinct(result, data)[i,], ".csv"))
}
