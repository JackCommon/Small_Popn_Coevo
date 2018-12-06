library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)
library(MuMIn)

#### Graphics functions & objects ####
# Next I've written some of simple functions that make the graph labels more 
# human-readable. There are also some objects that enable these functions to 
# work properly. *Importantly*, there are different ones for the different experiments, 
# which have been indicated in the comments.

# Lists of the original names (as found in the source data) of the levels of each factor. 
# These are used for the breaks() function in ggplots

## Timepoints are shared across experiments
OG_timepoints = c('t0', 't1', 't2', 't3', 't4', 't5')

timepoint_names_facet = list(
  't0' = '0 d.p.i.',
  't1' = '1 d.p.i.',
  't2' = '2 d.p.i.',
  't3' = '3 d.p.i.',
  't4' = '4 d.p.i.',
  't5' = '5 d.p.i.'
)

timepoint_names_legend = c(seq(0,20,1))

timepoint_labeller = function(variable, value) {
  return(timepoint_names_facet[value])
}

# These vectors just give tidy labels for the figures for experiment 3
oneclone_IDs <- c('M.1', 'M.2', 'M.3', 'M.4', 'M.5', 'M.6')
fiveclone_IDs <- c('5.1', '5.2', '5.3', '5.4', '5.5', '5.6')

plate_replicate_names_legend = c('1', '2', '3', '4', '5', '6')

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

pd = position_dodge(0.2)


##### Data ####
phage <- read.csv("./populations/original_data/phage_counts.csv", header = T)
phage <- select(phage, -cfu)
phage$ID %<>% as.factor()
phage %<>% na.exclude
phage$log.pfu <- log10(phage$pfu+1)

phage$timepoint%<>% relevel(ref="t16")
phage$timepoint%<>% relevel(ref="t15")
phage$timepoint%<>% relevel(ref="t14")
phage$timepoint%<>% relevel(ref="t13")
phage$timepoint%<>% relevel(ref="t12")
phage$timepoint %<>% relevel(ref="t11")
phage$timepoint %<>% relevel(ref="t10")
phage$timepoint %<>% relevel(ref="t9")
phage$timepoint %<>% relevel(ref="t8")
phage$timepoint %<>% relevel(ref="t7")
phage$timepoint %<>% relevel(ref="t6")
phage$timepoint %<>% relevel(ref="t5")
phage$timepoint %<>% relevel(ref="t4")
phage$timepoint %<>% relevel(ref="t3")
phage$timepoint %<>% relevel(ref="t2")
phage$timepoint %<>% relevel(ref="t1")
phage$timepoint %<>% relevel(ref="t0")

#### Subset dataframes by replicate
one <- filter(phage, ID=="1")
two <- filter(phage, ID=="2")
three <- filter(phage, ID=="3")
four <- filter(phage, ID=="4")
five <- filter(phage, ID=="5")
six <- filter(phage, ID=="6")
seven <- filter(phage, ID=="7")
eight <- filter(phage, ID=="8")
nine <- filter(phage, ID=="9")
ten <- filter(phage, ID=="10")
eleven <- filter(phage, ID=="11")
twelve <- filter(phage, ID=="12")

### Figures subsetted by replicate ID

#### One ####
one_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                         data=subset(one, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 1')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

one_mono_phage_plot

#### Two ####
two_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                             data=subset(two, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 2')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

two_mono_phage_plot

#### Three ####
three_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                             data=subset(three, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 3')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

three_mono_phage_plot

#### Four ####
four_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                             data=subset(four, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 4')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

four_mono_phage_plot

#### Five ####
five_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                             data=subset(five, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 5')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

five_mono_phage_plot

#### Six ####
six_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(six, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 6')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

six_mono_phage_plot

#### Seven ####
seven_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(seven, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 7')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

seven_mono_phage_plot

#### Eight ####
eight_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(eight, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 8')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

eight_mono_phage_plot

#### Nine ####
nine_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(nine, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 9')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

nine_mono_phage_plot

#### Ten ####
ten_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(ten, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 10')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

ten_mono_phage_plot

#### Eleven ####
eleven_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(eleven, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 11')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

eleven_mono_phage_plot

#### Twelve ####
twelve_mono_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(twelve, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicate 12')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12'),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

twelve_mono_phage_plot

grid <- plot_grid(one_mono_phage_plot, two_mono_phage_plot, three_mono_phage_plot,
          four_mono_phage_plot, five_mono_phage_plot, six_mono_phage_plot,
          seven_mono_phage_plot, eight_mono_phage_plot, nine_mono_phage_plot,
          ten_mono_phage_plot, eleven_mono_phage_plot, twelve_mono_phage_plot,
          nrow = 4, ncol=3)
quartz(grid)

#### Some more subset graphs ####
mono_t3 <- filter(phage, ID%in%c("2", "4", "6", "12"))

mono_t3_plot <- ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                      data=subset(mono_t3, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position = pd)+
  geom_path(stat='identity', position = pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicates extinct by 3 d.p.i.')+
  
  theme_cowplot()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12', "t13", "t14", "t15",
                            "t16"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15",
                            "16"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))
mono_t3_plot

mono_rest <- filter(phage, !(ID%in%c("2", "4", "6", "12")))

mono_rest_plot <- ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                       data=subset(mono_rest, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position = pd)+
  geom_path(stat='identity', position = pd, size=.5)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: replicates not extinct by 3 d.p.i.')+
  
  theme_cowplot()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12', "t13", "t14", "t15",
                            "t16"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15",
                            "16"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))
mono_rest_plot

mono_troopers <- filter(phage, (ID%in%c("10", "11")))

mono_troopers_plot <- ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                         data=subset(mono_troopers, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position = pd)+
  geom_path(stat='identity', position = pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone: absolute legends')+
  
  theme_cowplot()+
  scale_colour_discrete(name='Replicate',
                        breaks = oneclone_IDs,
                        labels = plate_replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5', 
                            't6', 't7', 't8', 't9', 't10',
                            't11', 't12', "t13", "t14", "t15",
                            "t16"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15",
                            "16"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))
mono_troopers_plot

ggsave("1-clone-t3-extinct.png", mono_t3_plot, path = "./presentation_figs/",
       device="png", dpi=300,
       width=28, height=15, units=c("cm"))

ggsave("1-clone-rest-noline.png", mono_rest_plot, path = "./presentation_figs/",
       device="png", dpi=300,
       width=28, height=15, units=c("cm"))

ggsave("1-clone-troopers.png", mono_troopers_plot, path = "./presentation_figs/",
       device="png", dpi=300,
       width=28, height=15, units=c("cm"))

mono_rest_plot <- mono_rest_plot + geom_vline(xintercept = c(2, 4, 6, 8),
                                              colour="red", linetype=2)

ggsave("1-clone-rest-WITHline.png", mono_rest_plot, path = "./presentation_figs/",
       device="png", dpi=300,
       width=28, height=15, units=c("cm"))
