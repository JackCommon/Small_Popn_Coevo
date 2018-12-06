#### Colonies_for_timeshift_gen.R ###
# Code to randomly generate 12 host isolates to select from the 24 original picks for each timepoint X replicate
library(magrittr)

selections <- data.frame(matrix(ncol=12, nrow=48))

for(i in seq(1,48)){
  selections[i,] <- sample(1:24, 12, replace=F)
}

selections

clip <- pipe('pbcopy', 'w')
write.table(selections, file=clip, row.names = F, col.names = F, sep="\t")
close(clip)
