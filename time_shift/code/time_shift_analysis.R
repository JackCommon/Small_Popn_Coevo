#### Small pop coevo - Time-shift analysis ####
# Created: 10/10/18 by Jack Common

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


#### Data Wrangling ####
# Original data is formatted as a matrix, which R will struggle with.
# Therefore, that data needs to be converted into an R-readable long-format

# Step 1: Melt each wide-formatted matrix into a long-format dataframe and joing
# them based on host timepoint

# 1-clone, replicate 1
one.pT1 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r1/1-clone_r1_phageT1.csv", header=T)
one.pT1 <- melt(one.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

one.pT3 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r1/1-clone_r1_phageT3.csv", header=T)
one.pT3 <- melt(one.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

one.pT5 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r1/1-clone_r1_phageT5.csv", header=T)
one.pT5 <- melt(one.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

one.pT7 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r1/1-clone_r1_phageT7.csv", header=T)
one.pT7 <- melt(one.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

one <- bind_rows(one.pT1, one.pT3, one.pT5, one.pT7)

# 1-clone, replicate 5
five.pT1 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r5/1-clone_r5_phageT1.csv", header=T)
five.pT1 <- melt(five.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

five.pT3 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r5/1-clone_r5_phageT3.csv", header=T)
five.pT3 <- melt(five.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

five.pT5 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r5/1-clone_r5_phageT5.csv", header=T)
five.pT5 <- melt(five.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

five.pT7 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r5/1-clone_r5_phageT7.csv", header=T)
five.pT7 <- melt(five.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

five <- bind_rows(five.pT1, five.pT3, five.pT5, five.pT7)

# 1-clone, replicate 7
seven.pT1 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r7/1-clone_r7_phageT1.csv", header=T)
seven.pT1 <- melt(seven.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

seven.pT3 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r7/1-clone_r7_phageT3.csv", header=T)
seven.pT3 <- melt(seven.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

seven.pT5 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r7/1-clone_r7_phageT5.csv", header=T)
seven.pT5 <- melt(seven.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

seven.pT7 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r7/1-clone_r7_phageT7.csv", header=T)
seven.pT7 <- melt(seven.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

seven <- bind_rows(seven.pT1, seven.pT3, seven.pT5, seven.pT7)

# 1-clone, replicate 8
eight.pT1 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r8/1-clone_r8_phageT1.csv", header=T)
eight.pT1 <- melt(eight.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eight.pT3 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r8/1-clone_r8_phageT3.csv", header=T)
eight.pT3 <- melt(eight.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eight.pT5 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r8/1-clone_r8_phageT5.csv", header=T)
eight.pT5 <- melt(eight.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eight.pT7 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r8/1-clone_r8_phageT7.csv", header=T)
eight.pT7 <- melt(eight.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eight <- bind_rows(eight.pT1, eight.pT3, eight.pT5, eight.pT7)

# 1-clone, replicate 9
nine.pT1 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r9/1-clone_r9_phageT1.csv", header=T)
nine.pT1 <- melt(nine.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nine.pT3 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r9/1-clone_r9_phageT3.csv", header=T)
nine.pT3 <- melt(nine.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nine.pT5 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r9/1-clone_r9_phageT5.csv", header=T)
nine.pT5 <- melt(nine.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nine.pT7 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r9/1-clone_r9_phageT7.csv", header=T)
nine.pT7 <- melt(nine.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nine <- bind_rows(nine.pT1, nine.pT3, nine.pT5, nine.pT7)

# 1-clone, replicate 10
ten.pT1 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r10/1-clone_r10_phageT1.csv", header=T)
ten.pT1 <- melt(ten.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

ten.pT3 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r10/1-clone_r10_phageT3.csv", header=T)
ten.pT3 <- melt(ten.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

ten.pT5 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r10/1-clone_r10_phageT5.csv", header=T)
ten.pT5 <- melt(ten.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

ten.pT7 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r10/1-clone_r10_phageT7.csv", header=T)
ten.pT7 <- melt(ten.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

ten <- bind_rows(ten.pT1, ten.pT3, ten.pT5, ten.pT7)

# 1-clone, replicate 11
eleven.pT1 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r11/1-clone_r11_phageT1.csv", header=T)
eleven.pT1 <- melt(eleven.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eleven.pT3 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r11/1-clone_r11_phageT3.csv", header=T)
eleven.pT3 <- melt(eleven.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eleven.pT5 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r11/1-clone_r11_phageT5.csv", header=T)
eleven.pT5 <- melt(eleven.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eleven.pT7 <- read.csv("./time_shift/original_data/infection_matrices/1-clone_r11/1-clone_r11_phageT7.csv", header=T)
eleven.pT7 <- melt(eleven.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

eleven <- bind_rows(eleven.pT1, eleven.pT3, eleven.pT5, eleven.pT7)

# Bind together each replicate dataframe for 1-clone treatment
monoclonal <- bind_rows(one, five, seven, eight, nine, ten, eleven)

# 5-clone, replicate 2
two.pT1 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r2/5-clone_r2_phageT1.csv", header=T)
two.pT1 <- melt(two.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

two.pT3 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r2/5-clone_r2_phageT3.csv", header=T)
two.pT3 <- melt(two.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

two.pT5 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r2/5-clone_r2_phageT5.csv", header=T)
two.pT5 <- melt(two.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

two.pT7 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r2/5-clone_r2_phageT7.csv", header=T)
two.pT7 <- melt(two.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

two <- bind_rows(two.pT1, two.pT3, two.pT5, two.pT7)

# 5-clone, replicate 3
three.pT1 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r3/5-clone_r3_phageT1.csv", header=T)
three.pT1 <- melt(three.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

three.pT3 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r3/5-clone_r3_phageT3.csv", header=T)
three.pT3 <- melt(three.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

three.pT5 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r3/5-clone_r3_phageT5.csv", header=T)
three.pT5 <- melt(three.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

three.pT7 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r3/5-clone_r3_phageT7.csv", header=T)
three.pT7 <- melt(three.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

three <- bind_rows(three.pT1, three.pT3, three.pT5, three.pT7)

# 5-clone, replicate 6
six.pT1 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r6/5-clone_r6_phageT1.csv", header=T)
six.pT1 <- melt(six.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

six.pT3 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r6/5-clone_r6_phageT3.csv", header=T)
six.pT3 <- melt(six.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

six.pT5 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r6/5-clone_r6_phageT5.csv", header=T)
six.pT5 <- melt(six.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

six.pT7 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r6/5-clone_r6_phageT7.csv", header=T)
six.pT7 <- melt(six.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

six <- bind_rows(six.pT1, six.pT3, six.pT5, six.pT7)

# 5-clone, replicate 7
sevenII.pT1 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r7/5-clone_r7_phageT1.csv", header=T)
sevenII.pT1 <- melt(sevenII.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

sevenII.pT3 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r7/5-clone_r7_phageT3.csv", header=T)
sevenII.pT3 <- melt(sevenII.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

sevenII.pT5 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r7/5-clone_r7_phageT5.csv", header=T)
sevenII.pT5 <- melt(sevenII.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

sevenII.pT7 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r7/5-clone_r7_phageT7.csv", header=T)
sevenII.pT7 <- melt(sevenII.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

sevenII <- bind_rows(sevenII.pT1, sevenII.pT3, sevenII.pT5, sevenII.pT7)

# 5-clone, replicate 9
nineII.pT1 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r9/5-clone_r9_phageT1.csv", header=T)
nineII.pT1 <- melt(nineII.pT1, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nineII.pT3 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r9/5-clone_r9_phageT3.csv", header=T)
nineII.pT3 <- melt(nineII.pT3, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nineII.pT5 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r9/5-clone_r9_phageT5.csv", header=T)
nineII.pT5 <- melt(nineII.pT5, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nineII.pT7 <- read.csv("./time_shift/original_data/infection_matrices/5-clone_r9/5-clone_r9_phageT7.csv", header=T)
nineII.pT7 <- melt(nineII.pT7, id.vars = c("Treatment", "Replicate", "Phage", "Phage.Background", "Host.Timepoint"))

nineII <- bind_rows(nineII.pT1, nineII.pT3, nineII.pT5, nineII.pT7)
 
# Bind together each replicate dataframe for 5-clone treatment
fiveclonal <- bind_rows(two, three, six, sevenII, nineII)

# Bind together each treatment dataframe
#data <- monoclonal
data <- bind_rows(monoclonal, fiveclonal)


# Rename some variables so they make more visual sense
data <- plyr::rename(data, c("variable"="Host.Genotype",
                             "value"="Infected",
                             "Phage"="Phage.Genotype"))
data$Replicate %<>% as.factor()
data$Phage.Genotype %<>% as.factor()
data$Treatment %<>% as.factor
#data %<>% filter(Treatment==c("1-clone")) #gets rid of weird blank cells 
                                          #include 5-clone for final data
data$Phage.Background %<>% as.factor()

# Add the inverse of the infectivity data to get resistance data
data$Resistant <- ifelse(data$Infected=="1",0,1)

Phage.Timepoint <- c(rep("1", 576),
                     rep("3", 576),
                     rep("5", 576),
                     rep("7", 576)) %>%
  rep(12)
data$Phage.Timepoint <- Phage.Timepoint %>%
  as.factor

Host.Genotype <- c()
for(i in seq(1,12)){
  Host.Genotype <- c(Host.Genotype, rep(i,48))
}

Host.Genotype %<>% rep(4) %>% rep(12)
data$Host.Genotype <- Host.Genotype %>% 
  as.factor

# Host.Background -> If phage are challenged against hosts from their past, present or future
# Phage.Background -> If hosts are challenged against phage from their past, present or future

# Finally, save the data as a CSV to the working directory
write.csv(file="./time_shift/original_data/infection_matrices/full_matrix.csv", data,
          row.names = F)

# Clear up intermediate dataframes
rm(monoclonal, fiveclonal,
   one, two, three, five, six, seven, sevenII, eight, nine, nineII, ten, eleven,
   one.pT1, one.pT3, one.pT5, one.pT7, two.pT1, two.pT3, two.pT5, two.pT7,
   three.pT1, three.pT3, three.pT5, three.pT7, five.pT1, five.pT3, five.pT5, five.pT7,
   six.pT1, six.pT3, six.pT5, six.pT7, seven.pT1, seven.pT3, seven.pT5, seven.pT7,
   sevenII.pT1, sevenII.pT3, sevenII.pT5, sevenII.pT7, eight.pT1, eight.pT3, eight.pT5, eight.pT7,
   nine.pT1, nine.pT3, nine.pT5, nine.pT7, nineII.pT1, nineII.pT3, nineII.pT5, nineII.pT7,
   ten.pT1, ten.pT3, ten.pT5, ten.pT7, eleven.pT1, eleven.pT3, eleven.pT5, eleven.pT7,
   Host.Genotype, Phage.Timepoint)


#### Reload the data ####
data <- read.csv("./time_shift/original_data/infection_matrices/full_matrix.csv")
data$Replicate %<>% as.factor()
data$Phage.Genotype %<>% as.factor()
data$Phage.Timepoint %<>% as.factor()
data$Host.Timepoint %<>% as.factor()
data$Host.Genotype %<>% as.factor()
#data$Infected %<>% as.numeric()
#data$SpacerNumber <- ifelse(is.na(data$SpacerNumber==T),0,data$SpacerNumber) #%>% as.factor
data %<>% na.exclude # include to get marginal and conditional R2 from models

monoclonal <- data %>% filter(Treatment==c("1-clone"))
fiveclonal <- data %>% filter(Treatment==c("5-clone"))

monoclonal2 <- filter(monoclonal, Host.Timepoint!="1") %>% 
  filter(Phage.Timepoint!="1")
fiveclonal2 <- filter(fiveclonal, Host.Timepoint!="1") %>% 
  filter(Phage.Timepoint!="1")

#### Analysis - 1-clone ####
# Host environment-only model
# Slope does not vary with respect to phage genotype
m1 <- glmer(Infected~Phage.Background+(1|Phage.Background),
            data=monoclonal,
            family=binomial())
summary(m1)
par(mfrow=c(2,2))
plot(m1)

# Genotype as a single random effect with no interaction

m2 <- glmer(Infected~Phage.Background+(1|Replicate), data=monoclonal, family=binomial)

summary(multcomp::glht(m2))
plot(m2)

anova(m1,m2, test="Chisq")
anova(m2, test="Chisq")
R2 <- r.squaredGLMM(m2)
R2
R2[1]/R2[3]*100

# Overall genotype x Environment model
# Slope varies for each phage genotype as a random effect

m3 <- glmer(Infected~Phage.Background+(1|Replicate)+(Phage.Timepoint|Host.Timepoint),
            data=monoclonal, family=binomial)
summary(m3)
anova(m3, test="Chisq")

plot(m3)

fixed3 <- fixef(m3)
logit2prob(fixed3[1])
logit2prob(fixed3[1]+fixed3[2])
logit2prob(fixed3[1]+fixed3[3])

anova(m1, m2, m3, test="Chisq")
logLik(m1)
logLik(m2)
logLik(m3)
AIC(m1, m2, m3) %>% compare_AICs()
# Although the heteroskedacity can't seem to shift, I'll move ahead with model 2 based on the anova,
# log-likelihood and AIC comparisons
summary(m3)
logit2prob(fixef(m3)[1])
logit2prob(fixef(m3)[1]+fixef(m3)[2])
logit2prob(fixef(m3)[1]+fixef(m3)[3])

multcomp::cftest(m3)
CIs.lmerMod <- confint.merMod(m3, parm="beta_")
CIs <- confint(m3, parm="beta_")
CIs.wald <- confint(m3, parm="beta_", method="Wald", level=0.95)
fixed3
CIs.wald

logit2prob(fixed3[1]+CIs.wald[2])
logit2prob(fixed3[1]+CIs.wald[5])

logit2prob(fixed3[1]+CIs.wald[3])
logit2prob(fixed3[1]+CIs.wald[6])

monoclonal$Phage.Background %<>% relevel(ref="Past")

#### Analysis - 5-clone ####
# Host environment-only model
# Slope does not vary with respect to phage genotype
m1 <- glmer(Infected~Phage.Background+(1|Phage.Background),
            data=fiveclonal,
            family=binomial())
summary(m1)
par(mfrow=c(2,2))
plot(m1)

# Genotype as a single random effect with no interaction

m2 <- glmer(Infected~Phage.Background+(1|Replicate), data=fiveclonal, family=binomial)

summary(multcomp::glht(m2))
plot(m2)

anova(m1,m2, test="Chisq")
anova(m2, test="Chisq")
R2 <- r.squaredGLMM(m2)
R2
R2[1]/R2[3]*100

# Overall genotype x Environment model
# Slope varies for each phage genotype as a random effect
m3 <- glmer(Infected~Phage.Background+(1|Replicate)+(Phage.Timepoint|Host.Timepoint),
            data=fiveclonal,
            family=binomial)
summary(multcomp::glht(m3))
anova(m3, test="Chisq")

plot(m3)

anova(m1, m2, m3, test="Chisq")
logLik(m1)
logLik(m2)
logLik(m3)
AIC(m1, m2, m3) %>% compare_AICs()

# Although the heteroskedacity can't seem to shift, I'll move ahead with model 2 based on the anova,
# log-likelihood and AIC comparisons
fixed3 <- fixef(m3)
logit2prob(fixed3[1])
logit2prob(fixed3[1]+fixed3[2])
logit2prob(fixed3[1]+fixed3[3])

CIs.wald <- confint(m3, parm="beta_", method="Wald", level=0.95)
fixed3
CIs.wald

logit2prob(fixed3[1]+CIs.wald[2])
logit2prob(fixed3[1]+CIs.wald[5])

logit2prob(fixed3[1]+CIs.wald[3])
logit2prob(fixed3[1]+CIs.wald[6])

fiveclonal$Phage.Background %<>% relevel(ref="Past")
#### Analysis - Comparing 

#### Analysis - Full data ####
m1 <- glmer(Infected~Phage.Background*Treatment+(1|Replicate)+(Phage.Timepoint|Host.Timepoint),
            data=data, family=binomial)
summary(m1)
anova(m1, test="Chisq")
r.squaredGLMM(m1)

plot(m1)

fixed1 <- fixef(m1)
logit2prob(fixed1[1])
logit2prob(fixed1[1]+fixed1[2])
logit2prob(fixed1[1]+fixed1[3])
logit2prob(fixed1[1]+fixed1[4])
logit2prob(fixed1[1]+fixed1[5])
logit2prob(fixed1[1]+fixed1[6])

multcomp::cftest(m1)
CIs.wald <- confint(m1, parm="beta_", method="Wald")
fixed1
CIs.wald

logit2prob(fixed1[1]+CIs.wald[2])
logit2prob(fixed1[1]+CIs.wald[8])

logit2prob(fixed1[1]+CIs.wald[3])
logit2prob(fixed1[1]+CIs.wald[9])

logit2prob(fixed1[1]+CIs.wald[4])
logit2prob(fixed1[1]+CIs.wald[10])

### Dropping T1 data

m2 <- glmer(Infected~Phage.Background*Treatment+(1|Replicate)+(Phage.Timepoint|Host.Timepoint),
            data=filter(data, Host.Timepoint!="1", Phage.Timepoint!="1"), family=binomial)
summary(m2)
anova(m2, test="Chisq")

plot(m2)

fixed2 <- fixef(m2)
logit2prob(fixed2[1])
logit2prob(fixed2[1]+fixed2[2])
logit2prob(fixed2[1]+fixed2[3])
logit2prob(fixed2[1]+fixed2[4])


multcomp::cftest(m2)
CIs.wald <- confint(m2, parm="beta_", method="Wald")
fixed2
CIs.wald

logit2prob(fixed2[1]+CIs.wald[2])
logit2prob(fixed2[1]+CIs.wald[8])

logit2prob(fixed2[1]+CIs.wald[3])
logit2prob(fixed2[1]+CIs.wald[9])

logit2prob(fixed2[1]+CIs.wald[4])
logit2prob(fixed2[1]+CIs.wald[10])

#### Figure - Timeshift means ####
timeshift_means <- read.csv("./time_shift/summary_data/timeshift_means_fullmod.csv", header=T)

timeshift_means$Environment %<>% relevel(ref="Future")
timeshift_means$Environment %<>% relevel(ref="Contemporary")
timeshift_means$Environment %<>% relevel(ref="Past")

p1 <- ggplot(aes(y=Mean.Infect, x=Environment, group=Group), data=timeshift_means)+
  geom_point(position = position_dodge(.5),
             size=3)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper), 
                position = position_dodge(.5),
                width=0, size=1)+
  geom_path(stat="identity", size=.8, linetype=2, position=position_dodge(.5))+
  facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Phage background", y="Proportion of\nhosts infected")+
  #ggtitle("Including data from 1 dpi")+
  theme(plot.title = element_text(colour = "red"),
        axis.title = element_text(face="bold", size=16),
        axis.text = element_text(size=14),
        legend.title = element_text(face='bold', size=14),
        legend.title.align = 0.5,
        legend.position = 'right',
        legend.key.width = unit(2, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size=14))+
  
  scale_y_continuous(breaks=c(seq(0, 1, 0.02)))+
  coord_cartesian(ylim=c(0,.1))+
  NULL
p1

ggsave("timeshift.png", p1, path="./figs/",
       dpi=300, device="png", width=25, height=12, unit=c("cm"))
#### Contemporary interactions ####

# data$Treatment %<>% relevel(ref="5-clone")
# data$Treatment %<>% relevel(ref="1-clone")

m3 <- glmer(Infected~Host.Timepoint*Treatment+(1|Replicate)+(1|Phage.Timepoint),
            data=filter(data, Phage.Background=="Contemporary"), family=binomial)
summary(m3)

fixed3 <- fixef(m3)
logit2prob(fixed3[1])
logit2prob(fixed3[1]+fixed3[2])
logit2prob(fixed3[1]+fixed3[3])
logit2prob(fixed3[1]+fixed3[4])

CIs.m3 <- confint(m3, method="Wald", parm="beta_")
CIs.m3

logit2prob(fixed3[1]+CIs.m3[2])
logit2prob(fixed3[1]+CIs.m3[10])

logit2prob(fixed3[1]+CIs.m3[3])
logit2prob(fixed3[1]+CIs.m3[11])

logit2prob(fixed3[1]+CIs.m3[4])
logit2prob(fixed3[1]+CIs.m3[12])

logit2prob(fixed3[1]+CIs.m3[5])
logit2prob(fixed3[1]+CIs.m3[13])



cont_means <- read.csv("./time_shift/summary_data/contemporary_means.csv", header=T)

p3 <- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=cont_means)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper),
                width=0, size=1)+
  geom_path(stat="identity", linetype=2, size=.8)+
  facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Proportion of\nhosts infected")+
  scale_x_discrete(breaks=c("t1", "t3", "t5", "t7"),
                   labels=c("1", "3", "5", "7"))+
  coord_cartesian(ylim=c(0,.5))+
  scale_y_continuous(breaks=c(seq(0,1,0.1)))+
  #scale_fill_discrete(name="Host timepoint\n(days post-infection)",
  #                   breaks=c("t1", "t4", "t9"),
  #                  labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  NULL

last_plot()

ggsave("contemporary_means.png", p3, path="./figs/",
       dpi=300, device="png", width=24, height=12, unit=c("cm"))
#### Figure - Contemporary means ####
cont_means <- read.csv("./time_shift/summary_data/contemporary_means.csv", header=T)

p2 <- ggplot(aes(y=Mean.Infect, x=Timepoint, group=Group), data=cont_means)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Infect.Lower, ymax=Infect.Upper),
                width=0, size=1)+
  geom_path(stat="identity", linetype=2, size=.8)+
  facet_wrap(~Treatment)+
  cowplot::theme_cowplot()+
  labs(x="Days post-infection", y="Proportion of\nhosts infected")+
  scale_x_discrete(breaks=c("t1", "t3", "t5", "t7"),
                   labels=c("1", "3", "5", "7"))+
  coord_cartesian(ylim=c(0,.5))+
  scale_y_continuous(breaks=c(seq(0,1,0.1)))+
  #scale_fill_discrete(name="Host timepoint\n(days post-infection)",
  #                   breaks=c("t1", "t4", "t9"),
  #                  labels=c("1", "4", "9"))+
  
  theme(axis.title = element_text(face="bold", size=16))+
  theme(axis.text = element_text(size=14))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(1, 'cm'))+
  theme(legend.text = element_text(size=14))+
  NULL

p2

ggsave("contemporary_means.png", p2, path="./figs/",
       dpi=300, device="png", width=25, height=12, unit=c("cm"))
