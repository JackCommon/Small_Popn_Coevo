### Small population coevolution experiment - phage titre analysis ####
# Created: 19/4/18 by Jack Common

rm(list=ls())

#### Dependencies ####
#install.packages("ggplot2")
#install.packages("scales")
#install.packages("reshape2")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("magrittr")
#install.packages("cowplot")
#install.packages("MuMIn")

library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)
library(MuMIn)

#### Stats functions ####
# Extracts the intercept coefficient (mean) and 95% CIs from the GLM objects, 
#with added functionality to copy those to the clipboard for easier input to summary dataframes

# Enable this if using Ubuntu
#clipboard <- function(x, sep="\t", row.names=FALSE, col.names=FALSE){
#    con <- pipe("xclip -selection clipboard -i", open="w")
#    write.table(x, con, sep=sep, row.names=row.names, col.names=col.names)
#    close(con)
#}

model_stats = function(model){
  sum = coef(model)
  conf = confint(model, level=c(0.95))
  print(c(sum[1]))
  print(c(conf[1,1]))
  print(c(conf[1,2]))
  
  stats = data.frame(conf[1,1], conf[1,2])
  #clipboard(stats)
  clip = pipe('pbcopy', 'w')
  write.table(stats, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Coefficients copied to the clipboard')
  
}

CopyCIs <- function(model){
  conf = logit2prob(confint(model, level=c(0.95)))
  print(c(conf[1,1]))
  print(c(conf[1,2]))
  CIs <- data.frame(conf[1,1], conf[1,2])
  
  #clipboard(CIs)
  clip = pipe('pbcopy', 'w')
  write.table(CIs, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Probability 95% CIs copied to the clipboard')
  
}

## Compare AIC values for model fit
# This function extracts the AIC for each GLM, and then compares the absolute relative 
# differences for each AIC to the model with the lowest AIC. This acts as a measure of 
# model fit. More can be found at: http://faculty.washington.edu/skalski/classes/QERM597/papers_xtra/Burnham%20and%20Anderson.pdf

compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

## Convert log-odds to probabilities from binomial GLM output
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

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

pd = position_dodge(0.1)

##### Data ####
phage <- read.csv("./phage_popns/original_data/phage_counts.csv", header = T)
phage <- select(phage, -cfu)
phage$ID %<>% as.factor()
phage %<>% na.exclude
phage$log.pfu <- log10(phage$pfu+1)

phage$timepoint %<>% relevel(ref="t15")
phage$timepoint %<>% relevel(ref="t14")
phage$timepoint %<>% relevel(ref="t13")
phage$timepoint %<>% relevel(ref="t12")
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

#### Raw phage titre by replicate plots ####
mono_phage_plot <- ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                         data=subset(phage, bottleneck == '1-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('1-clone')+
  
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
                            't11', 't12', "t13", "t14", "t15"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15"))+
  
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0, 12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

mono_phage_plot

fiveclone_phage_plot = ggplot(aes(y=log.pfu, x=timepoint, group=ID), 
                              data=subset(phage, bottleneck == '5-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("Log"*{}[bold("10")]*" p.f.u ml"*{}^{-1}*"")))+
  ggtitle('5-clone')+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = fiveclone_IDs,
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
                            't11', 't12', "t13", "t14", "t15"),
                   labels=c('0', '1', '2', '3', '4', '5',
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15"))+
  scale_y_continuous(breaks=c(seq(0,12,1)))+
  coord_cartesian(ylim=c(0,12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

fiveclone_phage_plot

## Arrange all the phage titer plots into one using plot_grid
mono_phage_plot <- mono_phage_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))
fiveclone_phage_plot <- fiveclone_phage_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))

sum.fig <- plot_grid(mono_phage_plot+labs(x='')+theme(legend.position = 'none'), 
                   fiveclone_phage_plot+theme(legend.position = 'none'),
                   ncol=1, nrow=2, align = "hv")
sum.fig

# Cowplot blocks ggsave so need to detach it
detach("package:cowplot")

ggsave("phage_fig.png", sum.fig, path="./figs/", device="png",
       width=28, height=20, unit=c("cm"), dpi=300)

#### Survival analysis #####

# Need to load these packages after the above as some important functions can be 
# blocked
#install.packages("survival")
#install.packages("rms")
#install.packages("car")
#install.packages("multcomp")
#install.packages("relaimpo")

library(survival)
library(rms)
library(car)
library(multcomp)
library(relaimpo)

phage<-read.csv("./phage_popns/original_data/survival_data.csv", header=T)
attach(phage)
names(phage)

#### Keplan-Meier graph ####
summary(KM<-survfit(Surv(time_to_death,status)~bottleneck))

jpeg("./figs/survplot.jpg", width=20, height=15, units="in", res=300)
par(mfrow=c(1,1), xpd=TRUE, oma=c(1.5,2.5,1,1), mai=c(1,1,1,1.2), bty="l", pty="s")

plot(survfit(Surv(phage$time_to_death,phage$status)~bottleneck), lty=c(1,3,5), lwd=c(5,5,5),
     ylab="", xlab="", axes=FALSE, ylim=c(0,1), xlim=c(0,16))

axis(1, tcl=-0.1, pos=0, cex.axis=1, lwd=c(3), cex.axis=2)
axis(1, at=8, lab="Days post-infection (d.p.i.)", tcl=0, line=2, cex.axis=3)

axis(2, tcl=-0.1, pos=-0, cex.axis=1, las=2, lwd=c(3), cex.axis = 2)
axis(2, at=0.5, lab="Proportion of phage\npopulations surviving", line=4, cex.axis=3, tcl=0)

legend(0.8,0.5, title=c("Bottleneck"),
       legend=c("1-clone", "5-clone"), 
       bty="o", lty=c(1,3,5), lwd=c(5,5,5), cex=3, adj=0)
dev.off()

#### Cox proportional hazards model ####
model3<-coxph(Surv(time_to_death,status)~bottleneck)
summary(model3)

model3$loglik

anova(model3)
tapply(predict(model3),bottleneck,mean)

