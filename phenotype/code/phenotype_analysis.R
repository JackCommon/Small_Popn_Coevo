### Small population covevolution ####
### Phenotype analysis ##
## Created: 24/5/18 by Jack Common

rm(list=ls())

### Dependencies ####
library(ggplot2)
library(scales)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)
library(lme4)

### Functions ####
# AIC comparison function
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

### Data ####
# Load in constituent datasets
t1 <- read.csv("./phenotype/original_data/phenotype_T1.csv")
t3 <- read.csv("./phenotype/original_data/phenotype_T3.csv")
t5 <- read.csv("./phenotype/original_data/phenotype_T5.csv")
t7 <- read.csv("./phenotype/original_data/phenotype_T7.csv")

#time.new <- c(rep(1,576), rep(3, 576), rep(5, 576), rep(7, 576))

# Bind to make the full dataset
data <- bind_rows(t1, t3, t5, t7)
data$Replicate %<>% as.factor()
data$Clone %<>% as.factor()
data$bottleneck %<>% as.factor()
#data$timepoint <- time.new
data$timepoint %<>% as.factor()
data <- melt(data, measure.vars = c('CRISPR', 'SM', 'Sensitive'))
data %<>% rename(phenotype=variable)

### Analysis ####
# Need a binomial glm(m)
# Make nested GLMMs with Replicate as a random factor in each
mod.null <- glmer(value~1+(1|Replicate), data=data,
                  family=binomial())
mod.1 <- glmer(value~phenotype+(1|Replicate), data=data,
                  family=binomial())
mod.2 <- glmer(value~bottleneck+(1|Replicate), data=data,
               family=binomial())
mod.3 <- glmer(value~timepoint+(1|Replicate), data=data,
               family=binomial())
mod.4 <- glmer(value~phenotype*bottleneck+(1|Replicate), data=data,
               family=binomial())
mod.5 <- glmer(value~phenotype*timepoint+(1|Replicate), data=data,
               family=binomial())
mod.6 <- glmer(value~bottleneck*timepoint+(1|Replicate), data=data,
               family=binomial())
mod.g <- glmer(value~phenotype*bottleneck*timepoint+(1|Replicate), data=data,
               family=binomial())
# Compare the fitted residuals among candidate models
plot(mod.null)
plot(mod.1)
plot(mod.2)
plot(mod.3)
plot(mod.4)
plot(mod.5)
plot(mod.6)
plot(mod.g)

# And check model fit with delta-AIC
AIC(mod.null, mod.1, mod.2, mod.3,
    mod.4, mod.5, mod.6, mod.g) %>% compare_AICs()

# Run a multi-model ANOVA to check for effect size
anova(mod.null, mod.1, mod.2, mod.3,
      mod.4, mod.5, mod.6, mod.g,
      test="Chisq")

# Looks like the global model is the best
# Least heteroskedacity, lowest AIC, and lowest log-likelihood
# Take a look at the model summary and coefficients
summary(mod.g)
# And just the fixed effects
fixed <- fixef(mod.g)
fixed %<>% as.data.frame()
fixed
logit2prob(fixed[1,1])
logit2prob(fixed[1,1]+fixed[4,1])
logit2prob(fixed[1,1]+fixed[5,1])
logit2prob(fixed[1,1]+fixed[6,1])
logit2prob(fixed[1,1]+fixed[7,1])

# Get 95% CIs (takes a while to calculate so find something else to do for 20 mins or so)
# parm is set as "beta_" so only the CIs of the fixed effects are calculated, which saves a little time
CI <- confint(mod.g, method="Wald", parm="beta_")
CI
logit2prob(CI[1]); logit2prob(CI[25])
logit2prob(CI[1]+CI[2]); logit2prob(CI[25]+CI[26]) 
logit2prob(CI[1]+CI[3]); logit2prob(CI[25]+CI[27])
logit2prob(CI[1]+CI[2]); logit2prob(CI[1]+CI[3])
logit2prob(CI[1]+CI[4]); logit2prob(CI[25]+CI[28])
logit2prob(CI[1]+CI[5]); logit2prob(CI[25]+CI[29])
logit2prob(CI[1]+CI[6]); logit2prob(CI[25]+CI[30])
logit2prob(CI[1]+CI[7]); logit2prob(CI[25]+CI[31])
length(CI)

test <- fixef(mod.g)
test %<>% as.data.frame()

data$bottleneck %<>% relevel(ref="1-clone")
data$phenotype %<>% relevel(ref="CRISPR")
data$timepoint %<>% relevel(ref="t5")
## Figures ####
pheno_sum <- read.csv('./phenotype/summary_data/phenotype_summary.csv')
pheno_sum %<>% na.exclude()
pheno_sum$variable %<>% relevel(., ref='Sensitive')
pheno_sum$variable %<>% relevel(., ref='SM')
pheno_sum$variable %<>% relevel(., ref='CRISPR')

pheno_sum_fig <- ggplot(aes(y=mean, x=bottleneck, group=variable), data=pheno_sum)+
  geom_bar(aes(fill=variable), stat='identity', size=3.5, position=position_dodge(1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0, size=0.7, position=position_dodge(1))+
  #geom_dl(aes(label=mean), method="last.points", position = position_dodge(1))+
  
  labs(x='Bottleneck', y='Proportion')+
  facet_grid(~timepoint)+
  
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5, size = 16))+
  theme(axis.title = element_text(face='bold', size=14, colour="black"),
        axis.text = element_text(colour="black"))+
  theme(strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(colour="black", fill="transparent"))+
  theme(legend.title = element_text(face='bold', size=12))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = "right")+
  theme(legend.key.width = unit(1, 'cm'))+
  theme(legend.key.height = unit(0.8, 'cm'))+
  theme(legend.text = element_text(size=11))+
  
  #scale_x_discrete(breaks=plate_OG_bottleneck, labels=plate_bottleneck_names_legend)+
  theme(axis.text = element_text(size=12))+
  
  scale_fill_manual(name='Genotype',
                    breaks = c('CRISPR', "Sensitive", "SM"),
                    labels = c("CRISPR", "Sensitive", "SM"),
                    values = c("#F8766D", "#619CFF", "#00BA38"))
pheno_sum_fig

ggsave("phenotype.png", pheno_sum_fig, path="./figs/",
       dpi=300, device="png",
       width=28, height=15, units=c("cm"))

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

p1 <- ggplot(aes(y=value, x=timepoint), data=data)+
  geom_jitter(aes(colour=phenotype), alpha=.2)+
  binomial_smooth(aes(colour=phenotype))+
  coord_cartesian(ylim = c(0,1))+
  facet_wrap(~bottleneck)+
  NULL
p1
