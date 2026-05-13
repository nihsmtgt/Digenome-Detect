
data=read_csv(paste0("./", output, ".csv"))  |> 
  arrange(desc(CLSCORE))  |> 
  filter(filtered==FALSE & pair_score<=lowestCL) 

FDR3 = tibble(Score = seq(500, 2000, 1)/100, Cleaved =0, nonCleaved =0)
for(i in 1:nrow(FDR3)){
  t = as.numeric(FDR3[i, "Score"] )
  FDR3[i, "Cleaved"] = nrow(filter(data, filtered == F & CLSCORE >= t & data=="RNP(+)"))
  FDR3[i, "nonCleaved"] = nrow(filter(data, filtered == F & CLSCORE >= t & data=="RNP(-)"))
}
FDR3 = mutate(FDR3, ratio = nonCleaved/Cleaved)
tl2 = min(FDR3[FDR3$ratio<=FDR,"Score"])
th2 = min(FDR3[FDR3$ratio==0,"Score"])

ggplot(FDR3)+
  geom_line(aes(Score, ratio), size = 0.5, color="#555555")+
  scale_x_continuous(limits = c(5,th2+5))+
  theme_classic(base_size = 12)+
  geom_vline(xintercept = tl2, linetype=2, col="#FF7777")+
  geom_hline(yintercept = FDR, linetype=2, col="#555555")+
  xlab("CL-score")+
  ylab("Non-cleaved / Cleaved")

write_csv(tibble(threshold_low=tl2, threshold_high = th2), paste0("./",output,"_threshold.csv"))
write_csv(FDR3, paste0("./",output,"_FDR.csv") )

