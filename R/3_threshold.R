
data=read_csv(paste0("./", output, ".csv"))  |> 
  arrange(desc(CLSCORE))  |> 
  filter(filtered==FALSE & pair_score<=7) 

data_f_NT=filter(data, data=="RNP(-)")
data_f_sample=filter(data, data=="RNP(+)")
data_f_NT=data_f_NT[order(data_f_NT$CLSCORE, decreasing = T), ]
detectF_score=distinct(data_f_sample, CLSCORE, .keep_all = T)

width=25
FDR2=foreach (i = seq(1,nrow(detectF_score)-width+1),.packages = "tidyverse", .combine = dplyr::bind_rows)%dopar% {
  temp=detectF_score[seq(i,i+width-1), ]
  temp2=filter(data_f_sample, CLSCORE<=temp[1,]$CLSCORE & CLSCORE>=temp[width,]$CLSCORE)
  temp3=filter(data_f_NT, CLSCORE<=temp[1,]$CLSCORE & CLSCORE>=temp[width,]$CLSCORE)
  temp4=tibble(RNP=nrow(temp2),NT=nrow(temp3), FDR=NT/RNP, CLSCOREmin=temp[width, ]$CLSCORE, CLSCOREmax=temp[1, ]$CLSCORE)
}
th2=max(data_f_NT$CLSCORE)
tl2=max(FDR2[FDR2$FDR>=1,"CLSCOREmin"])

write_csv(tibble(threshold_high=th2, threshold_low=tl2), paste0("./",output,"_threshold.csv"))
write_csv(FDR2, paste0("./",output,"_FDR.csv") )

ggplot(FDR2)+
  geom_point(aes(CLSCOREmin, FDR), size = 0.5, color="#555555")+
  scale_x_continuous(limits = c(5,th2+5))+
  theme_classic(base_size = 12)+
  xlab("CL-score")+
  ylab("Non-cleaved / Cleaved")

