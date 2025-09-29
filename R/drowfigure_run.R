library(eulerr)
library(Biostrings)
library(patchwork)
library(svglite)
library(tidyverse)
library(ggplot2)

output="./figure/Sample_DGS_SIDT1.70x_distinct_filt_toolkit"
if(!dir.exists("./figure")){
  dir.create("./figure")
}

#data
data=read_csv("./Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt_compared_seq_gg_align2.csv") |> 
  mutate(distance=if_else(distance==-1, ">6", as.character(distance))) |> 
  filter(TFtoolkit==T & data=="RNP(+)")
data_all=read_csv("./Sample_DGS_HAP_SIDT1.70x_toolkit2.5_toolkit_distinct.csv") 

gg=read_csv("D:/analysis/GGGenome/SIDT1_d6_new.csv")

#function
color="#444444"
ggbase=function(x) {ggplot(x)+
      theme_classic(base_size = 11, base_family = "sans", base_line_size = 0.3)+
      scale_color_manual(values = c("TRUE"="#777777", "FALSE"="#bdbdbd"))+
      scale_fill_manual(values = c("TRUE"="#777777", "FALSE"="#bdbdbd"))+
      guides(fill="none", color="none")+  
      theme(plot.tag = element_text(size=16), )
}

#ratio (distance<=6)
g=group_by(gg, PAMmatch) |> 
  count(distance) |> 
  replace_na(list("FALSE"=0, "TRUE"=0))

r=group_by(filter(data), PAMmatch) |> 
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

#distance histgram
figF=ggbase(filter(data))+
  geom_bar(aes(distance_2,fill=PAMmatch_2))+
  lims(x=factor(c(1:13)))+
  scale_y_continuous(expand = c(0,0), limits = c(0,120))+
  xlab("Distance")+
  ylab("Count")

figF

nrow(filter(data, distance_2>=7 & PAMmatch_2==F))/nrow(filter(data, distance_2>=7))

result=read_csv("./randomseqalign.csv")

D=ggbase(result)+
  geom_bar(aes(factor(distance_2, levels = c(0:20)), fill=PAMmatch_2))+
  xlab("Distance")+
  lims(x=factor(c(0:13)))+
  scale_y_continuous(expand = c(0,0), limits = c(0,5000))+
  ylab("Count")

D
nrow(filter(result, distance_2>=7))/nrow(result)*100

#figure2
Fig2=(plot_spacer()|E)/(figF|D)
Fig2
ggsave("./figure/Figure2.svg", width = 9*2, height = 9/1.618*2, units = "cm")


#scatter plot of sample
C=ggbase(mutate(filter(data_all, pair_score<2.5), data=if_else(data=="RNP(+)", "C", "NC")))+
  geom_jitter(aes(score, data), size=0.1, color=color)+
  scale_x_continuous(limits = c(0,20))+
  scale_y_discrete(expand = c(0,0.5))+
  guides(color="none")+
  ylab(NULL)+
  xlab("Cleavage score")+
  geom_vline(xintercept = 2.5, linetype=2, col="#FF7777")+
  geom_vline(xintercept = 1.5, linetype=2, col="#7777FF")
C
ggsave(paste0(output, "_scatterplot_toolkit.svg"), width=9, height = 9/1.618, units = "cm")


