#### paper_figs.R
#### Code to generate figures for paper on the small population coevolution experiment
#### Created: 11/1/19 by Jack Common

rm(list=ls())

#### Dependencies ####
library(tidyverse)
library(scales)
library(magrittr)
library(cowplot)
library(wesanderson)

#### Graphics functions ####
plate_OG_bottleneck = c('1-clone', '5-clone')

plate_bottleneck_names_facet = list(
  '1-clone'      = expression('1-clone'),
  '5-clone'          = expression('5-clone')
)

plate_bottleneck_names_legend = c(
  expression('1-clone'), 
  expression('5-clone')
)

# Function for experiment 3 bottleneck labelling
plate_bottleneck_labeller = function(variable, value) {
  return(plate_bottleneck_names_facet[value])
}

pal <- wes_palette("Royal1",4, type = "discrete")

#### Figure 1: Population dynamics ####
dynamics <- read.csv("./populations/original_data/phage_counts.csv", header = T)
dynamics$ID %<>% as.factor()

dynamics$timepoint %<>% relevel(ref="t16")
dynamics$timepoint %<>% relevel(ref="t15")
dynamics$timepoint %<>% relevel(ref="t14")
dynamics$timepoint %<>% relevel(ref="t13")
dynamics$timepoint %<>% relevel(ref="t12")
dynamics$timepoint %<>% relevel(ref="t11")
dynamics$timepoint %<>% relevel(ref="t10")
dynamics$timepoint %<>% relevel(ref="t9")
dynamics$timepoint %<>% relevel(ref="t8")
dynamics$timepoint %<>% relevel(ref="t7")
dynamics$timepoint %<>% relevel(ref="t6")
dynamics$timepoint %<>% relevel(ref="t5")
dynamics$timepoint %<>% relevel(ref="t4")
dynamics$timepoint %<>% relevel(ref="t3")
dynamics$timepoint %<>% relevel(ref="t2")
dynamics$timepoint %<>% relevel(ref="t1")
dynamics$timepoint %<>% relevel(ref="t0")

## Figure 1A: 1-clone dynamics

Fig1A <- ggplot(aes(y=pfu+1, x=timepoint, group=ID), data=subset(dynamics, bottleneck == '1-clone'))+
  geom_path(stat='identity')+
  geom_line(aes(y=cfu+1), colour="blue", linetype=1)+
  geom_hline(yintercept=1e+02, linetype=2, colour="black")+
  
  labs(x='Days post-infection', y=expression(bold("Pfu/Cfu ml"*{}^{-1}*"")))+
  ggtitle("1-clone")+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12', "t13", "t14", "t15",
                            "t16", "t17", "t18", "t19", "t20"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15",
                            "16", "17", "18", "19", "20"))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  NULL

#Fig1A

### Figure 1B: 5-clone dynamics

Fig1B <- ggplot(aes(y=pfu+1, x=timepoint, group=ID), data=subset(dynamics, bottleneck == '5-clone'))+
  geom_path(stat='identity')+
  geom_line(aes(y=cfu+1), colour="blue", linetype=1)+
  geom_hline(yintercept=1e+02, linetype=2, colour="black")+
  
  labs(x='Days post-infection', y=expression(bold("Pfu/Cfu ml"*{}^{-1}*"")))+
  ggtitle('5-clone')+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12', "t13", "t14", "t15",
                            "t16", "t17", "t18", "t19", "t20"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15",
                            "16", "17", "18", "19", "20"))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  NULL


#Fig1B

## Arrange Figure 1
Fig1A <- Fig1A + theme(plot.margin = unit(c(2,2,0,1), 'pt'))
Fig1B <- Fig1B + theme(plot.margin = unit(c(2,2,0,1), 'pt'))

Figure1 <- plot_grid(Fig1A+labs(x=""), Fig1B,
                     ncol=1, nrow=2, align = "hv", labels = c("A", "B"))
# Figure1

ggsave("Figure1.png", Figure1, path="./figs/paper_figs/", device="png",
       width=28, height=20, unit=c("cm"), dpi=300)


#### Figure 2: 1-clone timeshift results ####
timeshift_means <- read.csv("./time_shift/summary_data/timeshift_means_fullmod.csv", header=T)

timeshift_means$Environment %<>% relevel(ref="Future")
timeshift_means$Environment %<>% relevel(ref="Contemporary")
timeshift_means$Environment %<>% relevel(ref="Past")

## Figure 2A: 1-clone timeshift means

Fig2A <- ggplot(aes(y=Mean.Infect, x=Environment, group=Group), data=filter(timeshift_means, Treatment=="1-clone"))+
  geom_point(position = position_dodge(.5), size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2, position=position_dodge(.5))+
  # facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Phage background", y="Proportion of\nhosts infected")+
  #ggtitle("Including data from 1 dpi")+
  theme(plot.title = element_text(colour = "red"),
        axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.02)))+
  coord_cartesian(ylim=c(0,.1))+
  NULL
# Fig2A

## Figure 2B: 1-clone contemporary means

cont_means <- read.csv("./time_shift/summary_data/contemporary_means.csv", header=T)

Fig2B <- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=filter(cont_means, Treatment=="1-clone"))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper),
                width=0, size=1)+
  geom_path(stat="identity", linetype=2, size=.8)+
  #facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="", y="Proportion of\nhosts infected")+
  scale_x_discrete(breaks=c("t1", "t3", "t5", "t7"),
                   labels=c("1", "3", "5", "7"))+
  coord_cartesian(ylim=c(0,.5))+
  scale_y_continuous(breaks=c(seq(0,1,0.1)))+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14))+
  NULL

# Fig2B

## Figure 2C: 1-clone infectivity evolution
evolution <- read.csv("./infectivity/summary_data/infectivity_resistance_means.csv")
evolution$Timepoint %<>% as.factor
evolution$Group %<>% as.factor

Fig2C<- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=filter(evolution, Treatment=="1-clone"))+
  geom_point(position = position_dodge(.5), size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  #facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Infectivity")+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14, face="bold", margin=margin(2,0,2,0,"mm")))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,.5))+
  NULL
# Fig2C

## Figure 2D: 1-clone resistance evolution

Fig2D <- ggplot(aes(y=Mean.Resist, x=Timepoint, group=Group), data=filter(evolution, Treatment=="1-clone"))+
  geom_point(position = position_dodge(.5), size=3)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  #facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Resistance")+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14, face="bold", margin=margin(2,0,2,0,"mm")))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,1))+
  NULL

# Fig2D

## Arrange and save Figure 2

Figure2 <- plot_grid(Fig2A, Fig2B, Fig2C, Fig2D,
                     align="hv", labels = c("A", "B", "C", "D"))
# Figure2

ggsave("Figure2.png", Figure2, path="./figs/paper_figs/", device="png",
       dpi=300, width=32, height=20, units=c("cm"))

#### Figure 3: 5-clone timeshift results ####
## Figure 3A: 5-clone timeshift means

Fig3A <- ggplot(aes(y=Mean.Infect, x=Environment, group=Group), data=filter(timeshift_means, Treatment=="5-clone"))+
  geom_point(position = position_dodge(.5), size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2, position=position_dodge(.5))+
  # facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Phage background", y="Proportion of\nhosts infected")+
  #ggtitle("Including data from 1 dpi")+
  theme(plot.title = element_text(colour = "red"),
        axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.02)))+
  coord_cartesian(ylim=c(0,.1))+
  NULL
# Fig3A

## Figure 3B: 5-clone contemporary means

Fig3B <- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=filter(cont_means, Treatment=="5-clone"))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper),
                width=0, size=1)+
  geom_path(stat="identity", linetype=2, size=.8)+
  #facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="", y="Proportion of\nhosts infected")+
  scale_x_discrete(breaks=c("t1", "t3", "t5", "t7"),
                   labels=c("1", "3", "5", "7"))+
  coord_cartesian(ylim=c(0,.5))+
  scale_y_continuous(breaks=c(seq(0,1,0.1)))+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14))+
  NULL

# Fig3B

## Figure 3C: 5-clone infectivity evolution

Fig3C<- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=filter(evolution, Treatment=="5-clone"))+
  geom_point(position = position_dodge(.5), size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  #facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Infectivity")+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14, face="bold", margin=margin(2,0,2,0,"mm")))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,.5))+
  NULL
# Fig3C

## Figure 3D: 5-clone resistance evolution

Fig3D <- ggplot(aes(y=Mean.Resist, x=Timepoint, group=Group), data=filter(evolution, Treatment=="5-clone"))+
  geom_point(position = position_dodge(.5), size=3)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  #facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Resistance")+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14, face="bold", margin=margin(2,0,2,0,"mm")))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,1))+
  NULL

# Fig3D

## Arrange and save Figure 3

Figure3 <- plot_grid(Fig3A, Fig3B, Fig3C, Fig3D,
                     align="hv", labels = c("A", "B", "C", "D"))
# Figure3

ggsave("Figure3.png", Figure3, path="./figs/paper_figs/", device="png",
       dpi=300, width=32, height=20, units=c("cm"))
#### Figure 4: Phenotype and spacer results ####
phenotype <- read.csv('./phenotype/summary_data/phenotype_summary.csv')
phenotype %<>% na.exclude()
phenotype %<>% rename(Phenotype=variable,
                      Treatment=bottleneck,
                      Timepoint=timepoint)
phenotype$Timepoint %<>% as.factor
phenotype$Phenotype %<>% relevel(., ref='Sensitive')
phenotype$Phenotype %<>% relevel(., ref='SM')
phenotype$Phenotype %<>% relevel(., ref='CRISPR')

## Figure 4A: 1-clone phenotypes

Fig4A <- ggplot(aes(y=mean, x=Timepoint, group=Phenotype), data=filter(phenotype, Treatment=="1-clone"))+
  geom_col(aes(fill=Phenotype), colour="black", position=position_dodge(.9))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0, size=0.7, position=position_dodge(.9))+
  labs(x='Days post-infection', y='Relative frequency')+
  ggtitle("1-clone")+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face='bold', size=14, colour="black"),
        axis.text = element_text(colour="black", size=12),
        strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(colour="black", fill="transparent"),
        legend.position = "none")+

  scale_fill_manual(name='Immune\nphenotype',
                    breaks = c('CRISPR', "Sensitive", "SM"),
                    labels = c("CRISPR", "Sensitive", "SM"),
                    values = pal)+
  coord_cartesian(ylim=c(0,1))+
  NULL

#Fig4A

## Figure 4B: 5-clone phenotypes

Fig4B <- ggplot(aes(y=mean, x=Timepoint, group=Phenotype), data=filter(phenotype, Treatment=="5-clone"))+
  geom_col(aes(fill=Phenotype), colour="black", position=position_dodge(.9))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0, size=0.7, position=position_dodge(.9))+
  labs(x='Days post-infection', y='Relative frequency')+
  ggtitle("5-clone")+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face='bold', size=14, colour="black"),
        axis.text = element_text(colour="black", size=12),
        strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(colour="black", fill="transparent"),
        legend.position = "right",
        legend.title = element_text(face='bold', size=12),
        legend.title.align = 0.5,
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.8, 'cm'),
        legend.text = element_text(size=11))+
  
  scale_fill_manual(name='Immune\nphenotype',
                    breaks = c('CRISPR', "Sensitive", "SM"),
                    labels = c("CRISPR", "Sensitive", "SM"),
                    values = pal)+
  coord_cartesian(ylim=c(0,1))+
  NULL

# Fig4B

## Figure 4C: 1-clone spacers per clone

spacers <- read.csv("./spacers/summary_data/total_spacers.csv")
spacers$Timepoint %<>% as.factor()

Fig4C <- ggplot(aes(x=Timepoint, y=TotalSpacers), data=filter(spacers, Treatment=="1-clone"))+
  geom_col(position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, size=.8, position = position_dodge(.9))+
  coord_cartesian(ylim=c(0,4))+
  theme_cowplot()+
  #theme_black()+
  labs(x="Days post-infection", y="Total number of spacers")+
  ggtitle("1-clone")+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        plot.title = element_text(face="bold", hjust=0, size = 16))+
  NULL

# Fig4C

## Figure 4D: 5-clone spacers per clone

Fig4D <- ggplot(aes(x=Timepoint, y=TotalSpacers), data=filter(spacers, Treatment=="5-clone"))+
  geom_col(position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, size=.8, position = position_dodge(.9))+
  coord_cartesian(ylim=c(0,4))+
  theme_cowplot()+
  #theme_black()+
  labs(x="Days post-infection", y="Total number of spacers")+
  ggtitle("5-clone")+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        plot.title = element_text(face="bold", hjust=0, size = 16))+
  NULL

# Fig4D

## Arrange and save Figure 4

Fig4_left <- plot_grid(Fig4A+xlab(""), Fig4C,
                       align = "hv", ncol=1, labels=c("A", "C"))
Fig4_right <- plot_grid(Fig4B+labs(y="", x=""), Fig4D+ylab(""),
                       align = "hv", axis="l", ncol=1, labels=c("B", "D"))

Figure4 <- plot_grid(Fig4_left, Fig4_right, ncol=2,
                     align="hv", axis="l")
# Figure4

ggsave("Figure4.png", Figure4, path="./figs/paper_figs/", device="png",
       dpi=300, width=32, heigh=20, units = c("cm"))


#### Figure S1: Control population dynamics ####
control <- read.csv("./populations/original_data/control_counts.csv", header=T)
control %<>% select(-bottleneck, -treatment, -raw, -dilution)
control %<>% rename(Replicate=ID, Timepoint=timepoint)
control$Timepoint %<>% relevel(ref="t5")
control$Timepoint %<>% relevel(ref="t4")
control$Timepoint %<>% relevel(ref="t3")
control$Timepoint %<>% relevel(ref="t2")
control$Timepoint %<>% relevel(ref="t1")
control$Timepoint %<>% relevel(ref="t0")

FigS1 <- ggplot(aes(y=pfu+1, x=Timepoint, group=Replicate), data=control)+
  geom_path(stat='identity')+
  geom_line(aes(y=cfu+1), colour="blue", linetype=1)+
  geom_hline(yintercept=1e+02, linetype=2, colour="black")+
  
  labs(x='Days post-infection', y=expression(bold("Pfu/Cfu ml"*{}^{-1}*"")))+
  
  theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5'),
                   labels=c('0', '1', '2', '3', '4', '5'))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  NULL

# FigS1
ggsave("FigureS1.png", FigS1, path="./figs/paper_figs/", device="png",
       dpi=300, width=15, height = 10, unit=c("cm"))

