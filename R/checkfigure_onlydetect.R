library(eulerr)
library(Biostrings)
library(patchwork)
library(svglite)
library(tidyverse)
library(ggplot2)

output="./figure/Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt_compared_seq_gg_align2"
if(!dir.exists("./figure")){
  dir.create("./figure")
}

#data
data=read_csv("./Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt_seq_gg_align2.csv") |> 
  mutate(distance=if_else(distance==-1, ">6", as.character(distance))) 
data_all=read_csv("./Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt.csv") 

data_threshold=read_csv("./Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt_threshold.csv")
th=data_threshold$threshold_low
th2=data_threshold$threshold_high
data_FDR=read_csv("Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt_FDR.csv")
th_upper = filter(data_FDR, CLSCOREmin==th)$CLSCOREmax

data=mutate(data, TFdetect=if_else(CLSCORE>th & filtered == F & data == "RNP(+)", T, F)) |> 
  filter(TFdetect == T)

data_all=mutate(data_all, TFdetect=if_else(CLSCORE>th & filtered == F, T, F))

gg=read_csv("./SIDT1_d6_new.csv.gz")

#function
color="#444444"
ggbase=function(x) {ggplot(x)+
      theme_classic(base_size = 11, base_family = "sans", base_line_size = 0.3)+
      scale_color_manual(values = c("TRUE"="#777777", "FALSE"="#bdbdbd"))+
      scale_fill_manual(values = c("TRUE"="#777777", "FALSE"="#bdbdbd"))+
      guides(fill="none", color="none")+  
      theme(plot.tag = element_text(size=16), )
}

#CLScore description
prob = tibble(Number=0:10) |> 
  mutate(Probability = dpois(Number, 10/150))
ggbase(prob)+
  geom_col(aes(Number, Probability))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_continuous(breaks = 0:10)+
  xlab("Number of accumulated forward head")
ggsave("./figure/poissonprob.svg", width = 8, height = 8/1.618, units = "cm")

#CLscore boxplot
A=ggbase(filter(data, TFdetect ==T & PAMmatch==T & distance>0 & data=="RNP(+)"))+
  geom_boxplot(aes(factor(distance), CLSCORE, fill=PAMmatch), outlier.size = 0.25)+
  guides(fill="none")+
  xlab("Distance")+
  ylab("CL-score")+
  lims(x=factor(c(1:6)), y=c(0,300))

B=ggbase(filter(data, TFdetect ==T & PAMmatch==F & distance>0 & data=="RNP(+)"))+
  geom_boxplot(aes(factor(distance), CLSCORE, fill=PAMmatch), outlier.size = 0.25)+
  guides(fill="none")+
  xlab("Distance")+
  ylab("CL-score")+
  lims(x=factor(c(1:6)), y=c(0,300))

A+B

#scatter plot of sample
C=ggbase(mutate(filter(data_all, filtered==F & CLSCORE>7), data=if_else(data=="RNP(+)", "Cleaved", "Non-\ncleaved")))+
  geom_jitter(aes(CLSCORE, data), size=0.1, color=color)+
  scale_x_continuous(limits = c(7,50))+
  scale_y_discrete(expand = c(0,0.5))+
  guides(color="none")+
  ylab(NULL)+
  xlab("CL-score")+
  geom_vline(xintercept = th, linetype=2, col="#FF7777")+
  geom_vline(xintercept = th2, linetype=2, col="#7777FF")
C

data_all_2 = mutate(data_all, chr = str_remove(chr, "chr")) |> 
  filter(filtered==F & pair_score<=7)
C2=ggbase(filter(data_all_2, data=="RNP(+)" & CLSCORE>8 ))+
  geom_jitter(aes(factor(chr, levels = c(1:22, "X")), CLSCORE), size=0.5, color=color)+
  scale_y_continuous(expand = c(0,10), limits = c(0,300))+
  scale_x_discrete(limits=c(1:22, "X"))+
  guides(color="none")+
  ylab("CL-score")+
  xlab("Chromosome")
C2_2=ggbase(filter(data_all_2,  data=="RNP(-)" & CLSCORE>8) )+
  geom_jitter(aes(factor(chr, levels = c(1:22, "X")), CLSCORE), size=0.5, color=color)+
  scale_y_continuous(expand = c(0,10), limits = c(0,300))+
  scale_x_discrete(limits=c(1:22, "X"))+
  guides(color="none")+
  ylab("CL-score")+
  xlab("Chromosome")
C2/C2_2

Fig5=(A|B)/(C2/C2_2)
Fig5
ggsave(paste0(output, "_figure5.svg"), width=9*2, height = 9/1.618*2, units = "cm")


#FDR plot
D=ggbase(data_FDR)+
  geom_point(aes(CLSCOREmin, FDR), size = 0.2, color=color)+
  scale_x_continuous(limits = c(5,20))+
  xlab("CL-score")+
  ylab("Non-cleaved / Cleaved")+
  geom_vline(xintercept = th, linetype=2, col="#FF7777")+
  geom_vline(xintercept = th2, linetype=2, col="#7777FF")+
  geom_hline(yintercept = 1, linetype=1,linewidth=0.1)
D

#ratio (distance<=6)
g=group_by(gg, PAMmatch) |> 
  count(distance) |> 
  replace_na(list("FALSE"=0, "TRUE"=0))

r=group_by(filter(data,TFdetect==T & data=="RNP(+)"), PAMmatch) |> 
  count(distance)  |> 
  mutate(distance=if_else(distance==">6", "-1", distance)) |> 
  mutate(distance=as.double(distance))|> 
  replace_na(list("FALSE"=0, "TRUE"=0))

join=left_join(g,r, by=c("distance", "PAMmatch"))  |> 
  mutate(Ratio=n.y/n.x*100)

E=ggbase(filter(join,PAMmatch==T))+
  geom_col(aes(distance, Ratio))+
  lims(x=factor(c(1,2,3,4,5,6)))+
  scale_y_continuous(expand = c(0,0))+
  ylab("% Detected")+
  xlab("Distance")
E

#distance histogram
figF=ggbase(filter(data, TFdetect==T & data=="RNP(+)"))+
  geom_bar(aes(distance_2,fill=PAMmatch_2))+
  lims(x=factor(c(1:max(data$distance_2))))+
  scale_y_continuous(expand = c(0,0), limits = c(0,300))+
  xlab("Distance")+
  ylab("Count")

figF
ggsave(paste0(output, "_distancedist.svg"), width = 9, height = 9/1.618, units = "cm")

r2=group_by(filter(data,TFdetect==T & data=="RNP(+)"), PAMmatch_2) |> 
  count(distance_2) 


nrow(filter(data, TFdetect==T))
nrow(filter(data, TFdetect==T & distance_2>=7 & PAMmatch_2==F & data=="RNP(+)"))/nrow(filter(data, TFdetect==T & distance_2>=7 & data=="RNP(+)"))
nrow(filter(data, TFdetect==T & distance_2>=7 & data=="RNP(+)"))/nrow(filter(data, TFdetect==T & distance_2>=0 & data=="RNP(+)"))


#output
fig6=(D|C)/(E|figF)+plot_annotation(tag_levels = "A")

fig6
ggsave(paste0(output, "_figure6.svg"), width = 18, height = 9/1.618*2, units = "cm")

