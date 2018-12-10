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
  #scale_y_continuous(breaks=seq(0,1,.25))+
  coord_cartesian(ylim=c(0,4))+
  facet_wrap(~Treatment)+
  theme_cowplot()+
  labs(x="Days post-infection", y="Total number of spacers")+
  #scale_x_discrete(breaks=c("t1", "t4", "t9"),
   #                labels=c("1", "4", "9"))+
 # scale_fill_manual(name="Number of\nspacers",
  #                  values=pal)+
  
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
