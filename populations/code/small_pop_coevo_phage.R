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
  
  stats = data.frame(conf[1,1], conf[1,2])
  #clipboard(stats)
  clip = pipe('pbcopy', 'w')
  print(c(conf[1,2]))
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
phage <- read.csv("./populations/original_data/phage_counts.csv", header = T)
#phage <- select(phage, -cfu)
phage$ID %<>% as.factor()
#phage %<>% na.exclude
#phage$log.pfu <- log10(phage$pfu+1)

phage$timepoint %<>% relevel(ref="t16")
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
mono_phage_plot <- ggplot(aes(y=pfu+1, x=timepoint, group=ID), 
                         data=subset(phage, bottleneck == '1-clone'))+
  
  #geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  geom_line(aes(y=cfu+1), colour="blue", linetype=1)+
  
  labs(x='Days post-infection', y=expression(bold("Pfu/Cfu ml"*{}^{-1}*"")))+
  ggtitle('1-clone')+
  
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
                            "t16", "t17", "t18", "t19", "t20"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15",
                            "16", "17", "18", "19", "20"))+
  
  # scale_y_continuous(breaks=c(seq(0,12,1)))+
  # coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  NULL

mono_phage_plot

fiveclone_phage_plot = ggplot(aes(y=pfu+1, x=timepoint, group=ID), 
                              data=subset(phage, bottleneck == '5-clone'))+
  
 # geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  geom_line(aes(y=cfu+1), colour="blue", linetype=1)+
  
  labs(x='Days post-infection', y=expression(bold("Pfu/Cfu ml"*{}^{-1}*"")))+
  ggtitle('5-clone')+
  
  theme_cowplot()+
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
                            't11', 't12', "t13", "t14", "t15",
                            "t16", "t17", "t18", "t19", "t20"),
                   labels=c('0', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', '10',
                            '11', "12", "13", "14", "15",
                            "16", "17", "18", "19", "20"))+
  #scale_y_continuous(breaks=c(seq(0,12,1)))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  #coord_cartesian(ylim=c(0,12))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  NULL

fiveclone_phage_plot

## Arrange all the phage titer plots into one using plot_grid
mono_phage_plot <- mono_phage_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))
fiveclone_phage_plot <- fiveclone_phage_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))

sum.fig <- plot_grid(mono_phage_plot+labs(x='')+theme(legend.position = 'none',
                                                      panel.grid = element_blank(),
                                                      panel.grid.minor = element_blank(), 
                                                      panel.grid.major = element_blank(),
                                                      panel.background = element_blank(),
                                                      plot.background = element_rect(fill = "transparent",colour = NA)),
                   fiveclone_phage_plot+theme(legend.position = 'none'),
                   ncol=1, nrow=2, align = "hv")
sum.fig

ggsave("phage_fig.png", sum.fig, path="./figs/", device="png",
       width=28, height=20, unit=c("cm"), dpi=300)

#### Models of phage + bacteria titres ####
phage$bottleneck %<>% relevel(ref="5-clone")
phage$log.cfu <- log(phage$cfu+1)
phage$log.pfu <- log(phage$pfu+1)

m1 <- glm(log.pfu~log.cfu, data=phage)
m2 <- glm(log.pfu~log.cfu*bottleneck, data=phage)

par(mfrow=c(2,2))
plot(m1)
plot(m2)

AIC(m1, m2) %>% compare_AICs()
anova(m1, m2, test="Chisq")

summary(m2)
confint(m2)

association_plot <- ggplot(aes(y=pfu+1, x=cfu+1), data=phage)+
  geom_smooth(method="lm", formula = y~x, se=T, fullrange=T)+
  geom_point()+
  #stat_smooth(method="lm", formula = y~1, se=T, linetype=2, colour="grey2")+ # This gives a null slope
  facet_wrap(~bottleneck)+
  theme_cowplot()+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  scale_x_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  coord_cartesian(ylim=c(1,1e+10))+
  labs(x=expression(bold("C.f.u. ml"*{}^{-1}*"")), y=expression(bold("P.f.u. ml"*{}^{-1}*"")))+
  theme(strip.text.x = element_text(face="bold"))+
  NULL

quartz()
association_plot

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

phage<-read.csv("./populations/original_data/survival_data.csv", header=T)
attach(phage)
names(phage)

#### Keplan-Meier graph ####
summary(KM<-survfit(Surv(time_to_death,status)~bottleneck))
print(KM, print.rmean=TRUE)

jpeg("./figs/survplot.jpg", width=20, height=15, units="in", res=300)
par(mfrow=c(1,1), xpd=TRUE, oma=c(1.5,2.5,1,1), mai=c(1,1,1,1.2), bty="l", pty="s")

plot(survfit(Surv(phage$time_to_death,phage$status)~bottleneck), lty=c(1,3,5), lwd=c(5,5,5),
     ylab="", xlab="", axes=FALSE, ylim=c(0,1), xlim=c(0,20))

axis(1, tcl=-0.1, pos=0, cex.axis=1, lwd=c(3), cex.axis=2)
axis(1, at=10, lab="Days post-infection", tcl=0, line=2, cex.axis=3)

axis(2, tcl=-0.1, pos=-0, cex.axis=1, las=2, lwd=c(3), cex.axis = 2)
axis(2, at=0.5, lab="Proportion of phage\npopulations surviving", line=4, cex.axis=3, tcl=0)

legend(11.5,0.7, title=c("Bottleneck"),
       legend=c("1-clone", "5-clone"), 
       bty="o", lty=c(1,3,5), lwd=c(5,5,5), cex=3, adj=0)
dev.off()

#### Cox proportional hazards model ####
model3<-coxph(Surv(time_to_death,status)~bottleneck)
model3
summary(model3)

model3$loglik

anova(model3)
tapply(predict(model3),bottleneck,mean)


#### Binomial survival analysis ####
#install.packages("ggfortify")
#install.packages("survminer")
library(ggfortify)
library(survminer)

data <- read.csv("./populations/original_data/binomial_survival_data.csv", header = TRUE)

m.null <- glm(cbind(Alive, Dead)~1,family=binomial, data=data)
m1 <- glm(cbind(Alive, Dead)~Time,family=binomial, data=data)
m2 <- glm(cbind(Alive, Dead)~Treatment,family=binomial, data=data)
m.global <- glm(cbind(Alive, Dead)~Treatment*Time,family=binomial, data=data)

par(mfrow=c(2,2))
plot(m.null)
plot(m1)
plot(m2)
plot(m.global)

AIC(m.null, m1, m2, m.global) %>% compare_AICs()

summary(m.global)
anova(m.null, m1, m2, m.global,test="Chisq")
summary(glht(m.global, linfct = mcp(Treatment = "Tukey")))

m.global$coef[1] 
m.global$coef[2] 
m.global$coef[3] 

confint(m.global)

range(data$Time)
n<-seq(0,20,length.out=1000)
length(n)
intercept_TR1<-m.global$coef[1]
intercept_TR2<-m.global$coef[1]+m.global$coef[2]
lower_TR1<-confint(m.global)[1]
lower_TR2<-confint(m.global)[1]+confint(m.global)[2]
upper_TR1<-confint(m.global)[5]
upper_TR2<-confint(m.global)[5]+confint(m.global)[6]
#slope_Time<-m.global$coef[3]
slope_TR1 <- m.global$coef[3]
slope_TR2 <- m.global$coef[3]+m.global$coef[4]
slope_lower<-confint(m.global)[3]
slope_upper<-confint(m.global)[6]

#fitted_TR1<-exp(intercept_TR1+slope_Time*n)/(1+exp(intercept_TR1+slope_Time*n))
#fitted_TR2<-exp(intercept_TR2+slope_Time*n)/(1+exp(intercept_TR2+slope_Time*n))
fitted_TR1<-exp(intercept_TR1+slope_TR1*n)/(1+exp(intercept_TR1+slope_TR1*n))
fitted_TR2<-exp(intercept_TR2+slope_TR2*n)/(1+exp(intercept_TR2+slope_TR2*n))

fitted_lower_TR1<-exp(lower_TR1+slope_lower*n)/(1+exp(lower_TR1+slope_lower*n))
fitted_upper_TR1<-exp(upper_TR1+slope_upper*n)/(1+exp(upper_TR1+slope_upper*n))
fitted_lower_TR2<-exp(lower_TR2+slope_lower*n)/(1+exp(lower_TR2+slope_lower*n))
fitted_upper_TR2<-exp(upper_TR2+slope_upper*n)/(1+exp(upper_TR2+slope_upper*n))

par(mfrow=c(1,1))
plot(data$Time, data$Alive/(data$Alive+data$Dead), xlim=c(0,20))
points(fitted_TR1~n, type="l")
points(fitted_TR2~n, type="l", col="red")
points(fitted_lower_TR1~n, type="l", col="grey")
points(fitted_upper_TR1~n, type="l", col="grey")
points(fitted_lower_TR2~n, type="l", col="grey")
points(upper_TR2~n, type="l", col="grey")

par(las=1)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(4, 4, 4, 4)) 

jpeg("./figs/surv_proportional.jpg", width=20, height=15, units="in", res=300)

plot(data$Time, data$Alive/(data$Alive+data$Dead),axes=FALSE, type="n",ylab = "",xlab = "", xlim=c(0,20), cex.axis=2.5, cex.lab=2.5)
axis(side = 1, lwd = 3, cex.axis=2.5, pos = c(0,1), tck= -.01)
axis(side = 2, lwd = 3, cex.axis=2.5, pos = c(-.15,1))
title(ylab="Mean proportion surviving",line=6, cex.lab=2.5, tcl=0)
title(xlab="Days post-infection (d.p.i.)",line=4, cex.lab=2.5, tcl=0)

points(fitted_TR2~n, type="l", lwd=4,lty=1, col="#718B2E")
points(fitted_TR1~n, type="l", lwd=4,lty=1, col= "#D39200")

legend(1, 0.3, c("5-clone"),pch = c(19),lty = c(1), col=c('#718B2E'),lwd=4, cex=2,bty='n')
legend(1, 0.4, c("1-clone"),pch = c(19),lty = c(1), col=c('#D39200'),lwd=4, cex=2,bty='n')

points(data$Time[data$Treatment=="1-clone"], data$Alive[data$Treatment=="1-clone"]/(data$Alive[data$Treatment=="1-clone"]
                                                                +data$Dead[data$Treatment=="1-clone"]), pch=19, cex=2,col=c('#D39200'))
points(data$Time[data$Treatment=="5-clone"], data$Alive[data$Treatment=="5-clone"]/(data$Alive[data$Treatment=="5-clone"]
                                                                    +data$Dead[data$Treatment=="5-clone"]), pch=19, cex=2,col=c('#718B2E'))
dev.off()
