#### Small pop coevo exp - infectivity analysis ####
# Created: 10/5/18 by Jack Common
rm(list=ls())

#### Dependencies ####
#install.packages("nlme")
#install.packages("lme4")
#install.packages("MuMIn")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("magrittr")

library(nlme)
library(lme4)
library(MuMIn)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(reshape2)

#### Functions ####
# A function to quickly convert logit coefficients from a binomial GLM(M) 
# into more intuitive probability values
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

copy_stats = function(model1, model2){
  phage_means <- model.tables(aov(model1), "means")
  host_means <- model.tables(aov(model2), "means")
  genotype <- c()
  phage_coefs <- c()
  for ( i in seq(1,12)){
    genotype[i] <- i
    phage_coefs[i] <- phage_means$tables$Phage.Genotype[i]
  }
  phageDF <- data.frame(phage_coefs)
  
  host_coefs <- c()
  for ( i in seq(1,12)){
    genotype[i] <- i
    host_coefs[i] <- host_means$tables$Host.Genotype[i]
  }
  hostDF <- data.frame(host_coefs)
  
  coefs <- bind_cols(phageDF, hostDF)
  
  #clipboard(stats)
  clip <- pipe('pbcopy', 'w')
  write.table(coefs, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Coefficients copied to the clipboard')
  
}

host_stats = function(model){
  means = model.tables(aov(model), "means")
  genotype = c()
  coefs = c()
  for ( i in seq(1,12)){
    genotype[i] <- i
    coefs[i] <- means$tables$Host.Genotype[i]
  }
  df <- data.frame(coefs)
  
  #clipboard(stats)
  clip = pipe('pbcopy', 'w')
  write.table(df, file=clip, sep='\t', row.names = F, col.names = F)
  close(clip)
  print('Host resistance coefficients copied to the clipboard')
  
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

#### Data ####
# Import the full infection matrix dataset
data <- read.csv("./time_shift/original_data/infection_matrices/full_matrix.csv", header=T)
data$Replicate %<>% as.factor
data$Phage.Genotype %<>% as.factor
data$Phage.Timepoint %<>% as.factor
data$Host.Genotype %<>% as.factor()
data$Host.Timepoint %<>% as.factor()
data %<>% na.exclude()

monoclonal <- data %>% 
  filter(Treatment=="1-clone")
fiveclonal <- data %>% 
  filter(Treatment=="5-clone")

#### 1-clone - Proportion of hosts infected by each phage genotype ####
# Make subset dataframes of t1, t3, t5, and t7
h1 <- filter(monoclonal, Host.Timepoint=="1")
h3 <- filter(monoclonal, Host.Timepoint=="3")
h5 <- filter(monoclonal, Host.Timepoint=="5")
h7 <- filter(monoclonal, Host.Timepoint=="7")

p1 <- filter(monoclonal, Phage.Timepoint=="1")
p3 <- filter(monoclonal, Phage.Timepoint=="3")
p5 <- filter(monoclonal, Phage.Timepoint=="5")
p7 <- filter(monoclonal, Phage.Timepoint=="7")

# Then run a GLM to measure proportion of hosts infected by each phage genotype
# Change both the data and subset argument for timepoint and replicate, respectively
m1 <- glm(Infected~Phage.Genotype,
          data=p7, subset=c(Replicate=="11"),
          family = binomial())
m2 <- glm(Resistant~Host.Genotype,
          data=h7, subset=c(Replicate=="11"),
          family = binomial())

#### 5-clone - Proportion of hosts infected by each phage genotype ####
# Make subset dataframes of t1, t3, t5, and t7
h1 <- filter(fiveclonal, Host.Timepoint=="1")
h3 <- filter(fiveclonal, Host.Timepoint=="3")
h5 <- filter(fiveclonal, Host.Timepoint=="5")
h7 <- filter(fiveclonal, Host.Timepoint=="7")

p1 <- filter(fiveclonal, Phage.Timepoint=="1")
p3 <- filter(fiveclonal, Phage.Timepoint=="3")
p5 <- filter(fiveclonal, Phage.Timepoint=="5")
p7 <- filter(fiveclonal, Phage.Timepoint=="7")

# Then run a GLM to measure proportion of hosts infected by each phage genotype
# Change both the data and subset argument for timepoint and replicate, respectively
m1 <- glm(Infected~Phage.Genotype,
          data=p7, subset=c(Replicate=="9"),
          family = binomial())
m2 <- glm(Resistant~Host.Genotype,
          data=h7, subset=c(Replicate=="9"),
          family = binomial())

# Copy the coefficients to the clipboard
copy_stats(m1, m2)

#### Load in replicate data parts and bind ####
one <- read.csv("./infectivity/original_data/1-clone_1.csv")
five <- read.csv("./infectivity/original_data/1-clone_5.csv")
seven <- read.csv("./infectivity/original_data/1-clone_7.csv")
eight <- read.csv("./infectivity/original_data/1-clone_8.csv")
nine <- read.csv("./infectivity/original_data/1-clone_9.csv")
ten <- read.csv("./infectivity/original_data/1-clone_10.csv")
eleven <- read.csv("./infectivity/original_data/1-clone_11.csv")
two <- read.csv("./infectivity/original_data/5-clone_2.csv")
three <- read.csv("./infectivity/original_data/5-clone_3.csv")
six <- read.csv("./infectivity/original_data/5-clone_6.csv")
sevenII <- read.csv("./infectivity/original_data/5-clone_7.csv")
nineII <- read.csv("./infectivity/original_data/5-clone_9.csv")


# Bind and organise
data <- bind_rows(one, five, seven, eight, nine, ten, eleven,
                  two, three, six, sevenII, nineII)

data$replicate %<>% as.factor()
data$timepoint %<>% as.factor
data$phage %<>% as.factor()
data$host %<>% as.factor

# Rename some variables so they make more visual sense
data <- plyr::rename(data, c("treatment"="Treatment",
                             "replicate"="Replicate", 
                             "timepoint"="Timepoint",
                             "phage"="Phage.Genotype",
                             "host" = "Host.Genotype",
                             "infectivity"="Infected",
                             "resistance"="Resisted"))


# Save the bound data 
write.csv(file="./infectivity/summary_data/full_infect_resist_bygenotype.csv", data,
          row.names = F)

#### Load in data again ####
data <- read.csv("./infectivity/summary_data/full_infect_resist_bygenotype.csv")
data$Replicate %<>% as.factor
data$Phage.Genotype %<>% as.factor
data$Host.Genotype %<>% as.factor
data$Timepoint %<>% as.factor

data %<>% na.exclude()

data$Infected <- ifelse(data$Infected<0, 0, data$Infected)
data$Resisted <- ifelse(data$Resisted<0, 0, data$Resisted)

monoclonal <- data %>% 
  filter(Treatment=="1-clone")
fiveclonal <- data %>% 
  filter(Treatment=="5-clone")

#### 1-clone Mean infectivity ####
m1 <- glm(Infected~Timepoint,
            data=monoclonal,
            family=binomial(link="logit"))

m2 <- glmer(Infected~Timepoint+(1|Timepoint),
            data=monoclonal,
            family=binomial(link="logit"))

m3 <- glmer(Infected~Timepoint+(1|Replicate),
            data=monoclonal,
            family=binomial(link="logit"))


# Check fitted residuals...
# Heavy clustering
plot(m1)
# Much less clustering but some fanning
plot(m2)
plot(m3)


# Compare AICs
AICs <- AIC(m1, m2) %>% compare_AICs()
# Model 1 has the lowest AIC but they're roughly the same

# ANOVA to compare effect contributions
anova(m1, m2, test="Chisq")
drop1(m2, test="Chisq")
# Suggests that model 2 is the best at explaining variation in the data

# On balance and despite the increased heteroskedacity in the second model,
# I'll use this one moving forwards due to the AIC and Chisq tests

# Get the model coefficients and confidence intervals
summary(multcomp::glht(m3))
fixed <- coef(m1)
fixed

# Coefs
logit2prob(fixed[1])
logit2prob(fixed[1]+fixed[2])
logit2prob(fixed[1]+fixed[3])
logit2prob(fixed[1]+fixed[4])

monoclonal$Timepoint %<>% relevel(ref="7")
m1 <- glm(Infected~Timepoint,
          data=monoclonal,
          family=binomial(link="logit"))
CI <- confint(m1)
CI

logit2prob(CI[1]); logit2prob(CI[5])
logit2prob(fixed[1]+CI[2])
logit2prob(fixed[1]+CI[6])
logit2prob(fixed[1]+CI[3])
logit2prob(fixed[1]+CI[7])
logit2prob(fixed[1]+CI[4])
logit2prob(fixed[1]+CI[8])

#### 1-clone mean resistance ####
# Nested models with Timepoint as a fixed effect. Using a negative binomial family
# to account for zero inflation
m1 <- glm(Resisted~Timepoint,
               data=monoclonal,
               family=binomial(link="logit"))

m2 <- glmer(Resisted~Timepoint+(1|Replicate),
               data=monoclonal,
               family=binomial(link="logit"))


# Check fitted residuals...
# Heavy clustering
plot(m1)
# Much less clustering but some fanning
plot(m2)

# Compare AICs
AICs <- AIC(m1, m2) %>% compare_AICs()
# Model 2, with replicate as a random effect, has the lowest AIC

# ANOVA to compare effect contributions
anova(m1, m2, test="Chisq")
drop1(m1, test="Chisq")
# Suggests that model 2 is the best at explaining variation in the data

# Model 2 is probably the best
# Heteroskedacity in Model 2 suggests it is most appropriate
# Get the model coefficients and confidence intervals
summary(multcomp::glht(m2))
fixed <- fixef(m2)
fixed
logit2prob(fixed[1])
logit2prob(fixed[1]+fixed[2])
logit2prob(fixed[1]+fixed[3])
logit2prob(fixed[1]+fixed[4])

monoclonal$Timepoint %<>% relevel(ref="7")
monoclonal$Timepoint %<>% relevel(ref="5")
monoclonal$Timepoint %<>% relevel(ref="3")
monoclonal$Timepoint %<>% relevel(ref="1")

m2 <- glmer(Resisted~Timepoint+(1|Replicate),
            data=monoclonal,
            family=binomial(link="logit"))
CI <- confint(m2, parm="beta_")

logit2prob(CI)

logit2prob(fixed[1]+CI[2])
logit2prob(fixed[1]+CI[6])
logit2prob(fixed[1]+CI[3])
logit2prob(fixed[1]+CI[7])
logit2prob(fixed[1]+CI[4])
logit2prob(fixed[1]+CI[8])

logit2prob(CIs[1])+logit2prob(CIs[2])
logit2prob(CIs[1])+logit2prob(CIs[3])

logit2prob(CIs[4])
logit2prob(CIs[4])+logit2prob(CIs[5])
logit2prob(CIs[4])+logit2prob(CIs[6])
#### 5-clone Mean infectivity ####
m1 <- glm(Infected~Timepoint,
          data=fiveclonal,
          family=binomial(link="logit"))

m2 <- glmer(Infected~Timepoint+(1|Replicate),
          data=fiveclonal,
          family=binomial(link="logit"))


# Check fitted residuals...
# Heavy clustering
plot(m1)
# Much less clustering but some fanning
plot(m2)

# Compare AICs
AICs <- AIC(m1, m2) %>% compare_AICs()
#Both models are identical

# ANOVA to compare effect contributions
anova(m1, m2, test="Chisq")
drop1(m2, test="Chisq")
# Suggests that model 2 is the best at explaining variation in the data

# On balance and despite the increased heteroskedacity in the second model,
# I'll use this one moving forwards due to the AIC and Chisq tests

# Get the model coefficients and confidence intervals
summary(multcomp::glht(m1))
fixed <- coef(m1)
logit2prob(fixed)

# Coefs
logit2prob(fixed[1])
logit2prob(fixed[1]+fixed[2])
logit2prob(fixed[1]+fixed[3])
logit2prob(fixed[1]+fixed[4])

fiveclonal$Timepoint %<>% relevel(ref="1")
m1 <- glm(Infected~Timepoint,
          data=fiveclonal,
          family=binomial(link="logit"))

CI <- confint(m1)
CI

logit2prob(CI[1])
logit2prob(CI[5])

logit2prob(fixed[1]+CIs[2])
logit2prob(fixed[1]+CIs[6])
logit2prob(fixed[1]+CIs[3])
logit2prob(fixed[1]+CIs[7])
logit2prob(fixed[1]+CIs[4])
logit2prob(fixed[1]+CIs[8])
#### 5-clone mean resistance ####
# Nested models with Timepoint as a fixed effect. Using a negative binomial family
# to account for zero inflation
m1 <- glm(Resisted~Timepoint,
          data=fiveclonal,
          family=binomial(link="identity"))

m2 <- glmer(Resisted~Timepoint+(1|Replicate),
            data=fiveclonal,
            family=binomial(link="logit"))

summary(m2)
# Check fitted residuals...
# Heavy clustering
plot(m1)
# Much less clustering but some fanning
plot(m2)

# Compare AICs
AICs <- AIC(m1, m2) %>% compare_AICs()
# Model 2, with replicate as a random effect, has the lowest AIC, but both are basically the same (∆AIC<2)

# ANOVA to compare effect contributions
anova(m1, m2, test="Chisq")
drop1(m1, test="Chisq")
# Suggests that model 2 is the best at explaining variation in the data

# Model 2 is probably the best
# Heteroskedacity in Model 2 suggests it is most appropriate
# Get the model coefficients and confidence intervals
summary(multcomp::glht(m2))
fixed <- fixef(m2)
fixed
logit2prob(fixed[1])
logit2prob(fixed[1]+fixed[2])
logit2prob(fixed[1]+fixed[3])
logit2prob(fixed[1]+fixed[4])

data$Timepoint %<>% relevel(ref="7")
data$Timepoint %<>% relevel(ref="5")
data$Timepoint %<>% relevel(ref="3")
data$Timepoint %<>% relevel(ref="1")

m2 <- glmer(Resisted~Timepoint+(1|Replicate),
          data=data,
          family=binomial(link="logit"))
CI <- confint(m2, parm="beta_")

logit2prob(CI)

logit2prob(fixed[1]+CI[2])
logit2prob(fixed[1]+CI[6])
logit2prob(fixed[1]+CI[3])
logit2prob(fixed[1]+CI[7])
logit2prob(fixed[1]+CI[4])
logit2prob(fixed[1]+CI[8])
#### Figures ####

means <- read.csv("./infectivity/summary_data/infectivity_resistance_means.csv")
means$Timepoint %<>% as.factor
means$Group %<>% as.factor

infect_plot <- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=means)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Infectivity")+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14, face="bold", margin=margin(2,0,2,0,"mm")))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,.5))+
  NULL
infect_plot

resist_plot <- ggplot(aes(y=Mean.Resist, x=Timepoint, group=Group), data=means)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Resist.Lower, ymax=Resist.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2)+
  facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Resistance")+
  theme(axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14, face="bold", margin=margin(2,0,2,0,"mm")))+
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  coord_cartesian(ylim=c(0,1))+
  NULL
resist_plot

library(cowplot)
infect_resist_evo <- plot_grid(infect_plot+labs(x=""), resist_plot,
                               nrow=2, align = "hv", labels = c("A", "B"))

infect_resist_evo

ggsave("infect_resist_evolution.png", infect_resist_evo, path="./figs/",
       device="png", dpi=300, width=30, height = 24, units=c("cm"))
  
