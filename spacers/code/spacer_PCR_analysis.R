#### Spacer PCR analysis ####
#### Created: 10/12/18 by Jack Common

rm(list=ls())

#### Dependencies ####
library(tidyverse)
library(magrittr)
library(cowplot)
library(lme4)

#### Functions ####
# A function to quickly convert logit coefficients from a binomial GLM(M) 
# into more intuitive probability values
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

## Compare AIC values for Model fit
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

library(gridExtra)

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

#### Data ####
monoclonal <- read.csv("./spacers/original_data/1-clone.csv")
fiveclonal <- read.csv("./spacers/original_data/5-clone.csv")
data <- bind_rows(monoclonal, fiveclonal)
data %<>% na.exclude
data$Treatment %<>% as.factor
data$Timepoint %<>% as.factor
data$Replicate %<>% as.factor
data$Clone %<>% as.factor

monoclonal <- data %>% 
  filter(Treatment=="1-clone")
fiveclonal <- data %>% 
  filter(Treatment=="5-clone")

#### Total spacer analysis - 1-clone ####
m1 <- lmer(TotalSpacers~Timepoint+(1|Timepoint), data=monoclonal)
m2 <- lmer(TotalSpacers~Timepoint+(1|Replicate), data=monoclonal)
m3 <- lmer(TotalSpacers~Timepoint+(Timepoint|Replicate), data=monoclonal)

plot(m1)
plot(m2)
plot(m3)

AIC(m1, m2, m3) %>% compare_AICs
anova(m1, m2, m3, test="Chisq")
anova(m2, m3, test="Chisq")

summary(multcomp::glht(m3))
fixed <- fixef(m3)
CI <- confint(m3, parm="beta_")
CI

fixed[1]                       # T1 mean
CI[1]; CI[5]                   # T1 CIs
fixed[1]+fixed[2]              # T3 mean
fixed[1]+CI[2]; fixed[1]+CI[6] # T3 CIs
fixed[1]+fixed[3]              # T5 mean
fixed[1]+CI[3]; fixed[1]+CI[7] # T5 CIs
fixed[1]+fixed[4]              # T7 mean
fixed[1]+CI[4]; fixed[1]+CI[8] # T7 CIs

#### Total spacer analysis - 5-clone ####
m1 <- lmer(TotalSpacers~Timepoint+(1|Timepoint), data=fiveclonal)
m2 <- lmer(TotalSpacers~Timepoint+(1|Replicate), data=fiveclonal)
m3 <- lmer(TotalSpacers~Timepoint+(Timepoint|Replicate), data=fiveclonal)

plot(m1)
plot(m2)
plot(m3)

AIC(m1, m2, m3) %>% compare_AICs
anova(m1, m2, m3, test="Chisq")
anova(m2, m3, test="Chisq")

summary(multcomp::glht(m3))
fixed <- fixef(m3)
CI <- confint(m3, parm="beta_")
CI

fixed[1]                       # T1 mean
CI[1]; CI[5]                   # T1 CIs
fixed[1]+fixed[2]              # T3 mean
fixed[1]+CI[2]; fixed[1]+CI[6] # T3 CIs
fixed[1]+fixed[3]              # T5 mean
fixed[1]+CI[3]; fixed[1]+CI[7] # T5 CIs
fixed[1]+fixed[4]              # T7 mean
fixed[1]+CI[4]; fixed[1]+CI[8] # T7 CIs
#### Total spacer boxplot ####

total_boxplot <- ggplot(aes(x=Timepoint, y=TotalSpacers), data=data)+
  geom_boxplot()+
  facet_wrap(~Treatment)+
  theme_cowplot()+
  NULL

last_plot()

#### Total spacer figure ####
total_spacer <- read.csv("./spacers/summary_data/total_spacers.csv")
total_spacer$Timepoint %<>% as.factor()

total_plot <- ggplot(aes(x=Timepoint, y=TotalSpacers), data=total_spacer)+
  geom_col(position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, size=.8, position = position_dodge(.9))+
  coord_cartesian(ylim=c(0,4))+
  facet_wrap(~Treatment)+
  theme_cowplot()+
  #theme_black()+
  labs(x="Days post-infection", y="Total number of spacers")+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(face='bold', size=14),
        legend.title.align = 0.5,
        legend.position = 'right',
        legend.direction = "horizontal",
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size=14))+
  NULL
last_plot()

ggsave("total_spacers.png", total_plot, path="./figs/", device="png",
       dpi=300, width=25, height=13, unit=c("cm"))

#### Spacer number - data wrangling ####
data2 <- data %>% 
  gather(key="SpacerNumber", 
         measurement=c(Zero, One, Two, Three, Four, FiveOrMore), 
         factor_key = T)

monoclonal2 <- data2 %>% 
  filter(Treatment=="1-clone")
fiveclonal2 <- data2 %>% 
  filter(Treatment=="5-clone")

#### Spacer number - 1-clone analysis ####
m1 <- glmer(value~SpacerNumber+(1|SpacerNumber), data=monoclonal2, family=binomial())
m2 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=monoclonal2, family=binomial())
m3 <- glmer(value~SpacerNumber*Timepoint+(Timepoint|Replicate), data=monoclonal2, family=binomial())

plot(m1)
plot(m2)
plot(m3)

AIC(m1, m2, m3) %>% compare_AICs
anova(m1, m2, m3, test="Chisq")
anova(m2, m3, test="Chisq")

summary(multcomp::glht(m2))

monoclonal2$Timepoint %<>% relevel(ref="7")
m7 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=monoclonal2, family=binomial("logit"))

monoclonal2$Timepoint %<>% relevel(ref="5")
m5 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=monoclonal2, family=binomial("logit"))

monoclonal2$Timepoint %<>% relevel(ref="3")
m3 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=monoclonal2, family=binomial())

monoclonal2$Timepoint %<>% relevel(ref="1")
m1 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=monoclonal2, family=binomial())


fixed7 <- fixef(m7)
fixed5 <- fixef(m5)
fixed3 <- fixef(m3)
fixed1 <- fixef(m1)
CI_7 <- confint(m7, parm="beta_", method="Wald")
CI_5 <- confint(m5, parm="beta_", method="Wald")
CI_3 <- confint(m3, parm="beta_", method="Wald")
CI_1 <- confint(m1, parm="beta_", method="Wald")

copy_CI_main <- function(effects, confs, llevel, ulevel){
  mean <- logit2prob(effects[1])
  lower <- logit2prob(confs[llevel])
  upper <- logit2prob(confs[ulevel])
  df <- data.frame(mean, lower, upper)
  
  #clipboard(CIs)
  clip = pipe('pbcopy', 'w')
  write.table(df, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print(df)
  cat('\nMean and 95% CIs copied to the clipboard!')
  
}
copy_CI_diff <- function(fixed, CI, llevel, ulevel){
  mean <- logit2prob(fixed[1]+fixed[llevel])
  lower <- logit2prob(fixed[1]+CI[llevel])
  upper <- logit2prob(fixed[1]+CI[ulevel])
  df <- data.frame(mean, lower, upper)
  
  #clipboard(CIs)
  clip = pipe('pbcopy', 'w')
  write.table(df, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print(df)
  cat('\nMean and 95% CIs copied to the clipboard!')
  
}

copy_CI_main(fixed7, CI_7, 1, 25)
copy_CI_diff(fixed7, CI_7, 2, 26)
copy_CI_diff(fixed7, CI_7, 3, 27)
copy_CI_diff(fixed7, CI_7, 4, 28)
copy_CI_diff(fixed7, CI_7, 5, 29)
copy_CI_diff(fixed7, CI_7, 6, 30)

copy_CI_main(fixed5, CI_5, 1, 25)
copy_CI_diff(fixed5, CI_5, 2, 26)
copy_CI_diff(fixed5, CI_5, 3, 27)
copy_CI_diff(fixed5, CI_5, 4, 28)
copy_CI_diff(fixed5, CI_5, 5, 29)
copy_CI_diff(fixed5, CI_5, 6, 30)

copy_CI_main(fixed3, CI_3, 1, 25)
copy_CI_diff(fixed3, CI_3, 2, 26)
copy_CI_diff(fixed3, CI_3, 3, 27)
copy_CI_diff(fixed3, CI_3, 4, 28)
copy_CI_diff(fixed3, CI_3, 5, 29)
copy_CI_diff(fixed3, CI_3, 6, 30)

copy_CI_main(fixed1, CI_1, 1, 25)
copy_CI_diff(fixed1, CI_1, 2, 26)
copy_CI_diff(fixed1, CI_1, 3, 27)
copy_CI_diff(fixed1, CI_1, 4, 28)
copy_CI_diff(fixed1, CI_1, 5, 29)
copy_CI_diff(fixed1, CI_1, 6, 30)

#### Spacer number - 5-clone analysis ####
m1 <- glmer(value~SpacerNumber+(1|SpacerNumber), data=fiveclonal2, family=binomial())
m2 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=fiveclonal2, family=binomial())
m3 <- glmer(value~SpacerNumber*Timepoint+(Timepoint|Replicate), data=fiveclonal2, family=binomial())

plot(m1)
plot(m2)
plot(m3)

AIC(m1, m2, m3) %>% compare_AICs
anova(m1, m2, m3, test="Chisq")
anova(m2, m3, test="Chisq")

summary(multcomp::glht(m2))

fiveclonal2$Timepoint %<>% relevel(ref="7")
m7 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=fiveclonal2, family=binomial("logit"))

fiveclonal2$Timepoint %<>% relevel(ref="5")
m5 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=fiveclonal2, family=binomial("logit"))

fiveclonal2$Timepoint %<>% relevel(ref="3")
m3 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=fiveclonal2, family=binomial())

fiveclonal2$Timepoint %<>% relevel(ref="1")
m1 <- glmer(value~SpacerNumber*Timepoint+(1|Replicate), data=fiveclonal2, family=binomial())


fixed7 <- fixef(m7)
fixed5 <- fixef(m5)
fixed3 <- fixef(m3)
fixed1 <- fixef(m1)
CI_7 <- confint(m7, parm="beta_", method="Wald")
CI_5 <- confint(m5, parm="beta_", method="Wald")
CI_3 <- confint(m3, parm="beta_", method="Wald")
CI_1 <- confint(m1, parm="beta_", method="Wald")

copy_CI_main <- function(effects, confs, llevel, ulevel){
  mean <- logit2prob(effects[1])
  lower <- logit2prob(confs[llevel])
  upper <- logit2prob(confs[ulevel])
  df <- data.frame(mean, lower, upper)
  
  #clipboard(CIs)
  clip = pipe('pbcopy', 'w')
  write.table(df, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print(df)
  cat('\nMean and 95% CIs copied to the clipboard!')
  
}
copy_CI_diff <- function(fixed, CI, llevel, ulevel){
  mean <- logit2prob(fixed[1]+fixed[llevel])
  lower <- logit2prob(fixed[1]+CI[llevel])
  upper <- logit2prob(fixed[1]+CI[ulevel])
  df <- data.frame(mean, lower, upper)
  
  #clipboard(CIs)
  clip = pipe('pbcopy', 'w')
  write.table(df, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print(df)
  cat('\nMean and 95% CIs copied to the clipboard!')
  
}

copy_CI_main(fixed7, CI_7, 1, 25)
copy_CI_diff(fixed7, CI_7, 2, 26)
copy_CI_diff(fixed7, CI_7, 3, 27)
copy_CI_diff(fixed7, CI_7, 4, 28)
copy_CI_diff(fixed7, CI_7, 5, 29)
copy_CI_diff(fixed7, CI_7, 6, 30)

copy_CI_main(fixed5, CI_5, 1, 25)
copy_CI_diff(fixed5, CI_5, 2, 26)
copy_CI_diff(fixed5, CI_5, 3, 27)
copy_CI_diff(fixed5, CI_5, 4, 28)
copy_CI_diff(fixed5, CI_5, 5, 29)
copy_CI_diff(fixed5, CI_5, 6, 30)

copy_CI_main(fixed3, CI_3, 1, 25)
copy_CI_diff(fixed3, CI_3, 2, 26)
copy_CI_diff(fixed3, CI_3, 3, 27)
copy_CI_diff(fixed3, CI_3, 4, 28)
copy_CI_diff(fixed3, CI_3, 5, 29)
copy_CI_diff(fixed3, CI_3, 6, 30)

copy_CI_main(fixed1, CI_1, 1, 25)
copy_CI_diff(fixed1, CI_1, 2, 26)
copy_CI_diff(fixed1, CI_1, 3, 27)
copy_CI_diff(fixed1, CI_1, 4, 28)
copy_CI_diff(fixed1, CI_1, 5, 29)
copy_CI_diff(fixed1, CI_1, 6, 30)

#### Spacer number - figures ####
spacer_numbers <- read.csv("./spacers/summary_data/spacer_numbers.csv")
spacer_numbers$Timepoint %<>% as.factor
spacer_numbers$SpacerNumber %<>% relevel(ref="FiveOrMore")
spacer_numbers$SpacerNumber %<>% relevel(ref="Four")
spacer_numbers$SpacerNumber %<>% relevel(ref="Three")
spacer_numbers$SpacerNumber %<>% relevel(ref="Two")
spacer_numbers$SpacerNumber %<>% relevel(ref="One")
spacer_numbers$SpacerNumber %<>% relevel(ref="Zero")

spacer_no <- ggplot(aes(x=Timepoint, y=Mean), data=spacer_numbers)+
  geom_col(position=position_dodge(.9), aes(fill=SpacerNumber), colour="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper, group=SpacerNumber), width=0, size=.8, position = position_dodge(.9))+
  coord_cartesian(ylim=c(0,1))+
  facet_wrap(~Treatment)+
  theme_cowplot()+
  #theme_black()+
  labs(x="Days post-infection", y="Relative frequency")+
  
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(face='bold', size=14),
        legend.title.align = 0.5,
        legend.position = 'none',
        legend.direction = "horizontal",
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size=14))+
  NULL
last_plot()
