#### Time-shift analysis ####
#### Small Population Coevolution Experiment ###
# Created: 15/5/18 by Jack Common

#### Dependencies ####
#install.packages("pwr")
#install.packages("compute.es")
library(pwr)        # for power analyses
library(compute.es) # to derive reference D

#### Power analysis ####
# Need to decide on what sample size will be sufficient to detect a given effect of coevolution
# Take a look at the ./time_shift/original_data/replicate_decision.xlsx file and the project log.
# Can either sample from T1, 3 and 5 with a N = 8 and 9 in 1-clone and 5-clone respectively
# OR sample from T1, 3, 5, and 7 with N = 7 and 5 in 1-clone and 5-clone respectively
# Need to do a power analysis to figure out if the lower N option is okay

# Some reference Cohen's D from:
# Morley et al 2018, representing a weak FSD effect:
es1 <- fes(f=1.3143, n.1=8, n.2=8)
es1$d
# And Hall et al 2011, representing a stronger FSD effect:
es2 <- fes(f=4.89, n.1=6, n.2=6)
es2$d

## Add the reference Cohen's D 
## Small ES
pwr.f2.test(u=1, v=7, f2 = 0.57)$power #1-clone, T5
pwr.f2.test(u=1, v=6, f2 = 0.57)$power #1-clone, T7
pwr.f2.test(u=1, v=8, f2 = 0.57)$power #5-clone, T5
pwr.f2.test(u=1, v=4, f2 = 0.57)$power #5-clone, T7

## Large ES
pwr.f2.test(u=1, v=7, f2 = 1.28)$power #1-clone, T5
pwr.f2.test(u=1, v=6, f2 = 1.28)$power #1-clone, T7
pwr.f2.test(u=1, v=8, f2 = 1.28)$power #5-clone, T5
pwr.f2.test(u=1, v=4, f2 = 1.28)$power #5-clone, T7

# Looks like the smallest N option (N=5 in 5-clone) will have sufficient power to detect
# even a small effect size
