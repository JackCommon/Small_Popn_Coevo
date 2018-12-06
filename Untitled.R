m3 <- glmer(Infected~Phage.Background+(1|Replicate)+(Phage.Timepoint|Host.Timepoint),
            data=fiveclonal2, family=binomial)
summary(m3)
anova(m3, test="Chisq")

plot(m3)

fixed <- fixef(m3)
logit2prob(fixed[1])
logit2prob(fixed[1]+fixed[2])
logit2prob(fixed[1]+fixed[3])
