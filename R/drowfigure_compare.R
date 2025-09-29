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

#compared
data_threshold=read_csv("./Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt_threshold.csv")
th=data_threshold$threshold_low

data=read_csv("./Sample_DGS_HAP_SIDT1.70x_toolkit2.5_distinct_filt_compared_seq_gg_align2.csv") |> 
  mutate(distance=if_else(distance==-1, ">6", as.character(distance))) |> 
  filter(data=="RNP(+)") |> 
  mutate(TFdetect=if_else(detect>th, T, F))

toolkit_only=function(x) filter(x, TFtoolkit==T & TFdetect==F)
detect_only=function(x) filter(x, TFtoolkit==F & TFdetect==T)
both=function(x) filter(x, TFtoolkit==T & TFdetect==T)
color="#444444"
ggbase=function(x) {ggplot(x)+
      theme_classic(base_size = 11, base_family = "sans", base_line_size = 0.3)+
      scale_color_manual(values = c("TRUE"="#636363", "FALSE"="#bdbdbd"))+
      scale_fill_manual(values = c("TRUE"="#636363", "FALSE"="#bdbdbd"))+
      guides(fill="none", color="none")
}

#eular diagram
venndata=select(data, TFtoolkit, TFdetect) 
euler(venndata, shape="ellipse") |> 
  plot(quantities=TRUE)

#distance_histgram
distbar=function(x) {ggbase(x)+
    geom_bar(aes(factor(distance_2, levels = c(0:20)), fill=PAMmatch_2))+
    scale_y_continuous(expand = c(0,0), limits = c(0,220))+
    lims(x=factor(c(1:13)))+
    guides(fill="none")+
    xlab("Distance")+
    ylab("Count")
}
bar_toolkit=distbar(toolkit_only(data))
bar_detect=distbar(detect_only(data))
bar_both=distbar(both(data))
A=bar_both+bar_toolkit+bar_detect+
  plot_layout(ncol=1)
A
ggsave(paste0(output, "_distbarall.svg"), width = 9, height = 9/1.618*3, units = "cm")

nrow(filter(toolkit_only(data), distance_2>=7 & PAMmatch_2==F))/nrow(filter(toolkit_only(data),distance_2>=7))

nrow(filter(detect_only(data), distance_2<7))/nrow(filter(detect_only(data)))
nrow(filter(detect_only(data), distance_2>=7 & PAMmatch_2==T))/nrow(filter(detect_only(data), distance_2>=7))



detect_d10=filter(detect_only(data), distance_2>=9)
toolkit_d10=filter(toolkit_only(data), distance_2>=9)
write_csv(detect_d10, "./test.csv")

#CFD score
data_CFD=mutate(data, analysis=if_else(TFtoolkit==T&TFdetect==T, "rd", if_else(TFtoolkit==T, "r", "d")))

ggbase(filter(data_CFD, CFD!=-1, distance_2>=9))+
  geom_boxplot(aes(analysis, CFD), position = position_dodge2(preserve = "single", reverse = T))+
  lims(x=factor(c("r","rd","d")))

ggbase(filter(data_CFD, CFD!=-1, distance_2>=9))+
  geom_jitter(aes(analysis, CFD))+
  lims(x=factor(c("r","rd","d")))

ggsave(paste0(output, "_CFDd7.svg"), width = 3, height = 1.618*3)  

a=filter(data, TFtoolkit==T & TFdetect==F & distance_2<=6)

#seed region distance
seeddistbox2=function(x) {ggbase(filter(x, distance==">6"))+
    geom_boxplot(aes(factor(distance_2, levels = c(0:20)), seeddist/distance_2, fill=factor(PAMmatch_2, levels = c(F,T))), position = position_dodge2(preserve = "single", reverse = T))+
    scale_y_continuous(limits = c(0,1), expand = c(0,0))+
    lims(x=factor(c(7:max(data$distance_2))))+
    xlab("Distance")+
    guides(fill="none")+
    geom_hline(yintercept = 0.6)
}

seedbox2_toolkit=seeddistbox2(toolkit_only(data))
seedbox2_detect=seeddistbox2(detect_only(data))
seedbox2_both=seeddistbox2(both(data))

C=seedbox2_both+seedbox2_toolkit+seedbox2_detect+
  plot_layout(ncol=1)
C
ggsave(paste0(output, "_seeddistbox.svg"), width = 3, height = 3/1.618*3)

#seed distance distribution
ggbase(toolkit_only(filter(data, distance_2>6)))+
  geom_bar(aes(seeddist, fill = PAMmatch_2))+
  xlim(c(0,8))

ggbase(detect_only(filter(data, distance_2>6)))+
  geom_bar(aes(seeddist, fill = PAMmatch_2))+
  xlim(c(0,8))

ggbase(detect_only(filter(data, distance_2>6 & PAMmatch_2==T)))+
  geom_point(aes(detect, seeddist))


