library(gridExtra); library(grid); require(janitor);
library(MuMIn);



#***********kinetics/decay
### decay ANOVA  ###
#Now read decay data and merge -- decay column called lambda

options(na.action = "na.fail")
#doesn't work inside lm function but can cut and paste manually
#after removing invariant sites, which are 160, 303, 338
paste(colnames(pos_sel_lam)[3:ncol(pos_sel_lam)-1], collapse=" + "   )
#IF sites change copy/paste from result of above command

glb1 <- lm(lambda ~ s41 + s93 + s102 + s160 + s177 + s189 + s261 + s291 + s389 + s477, data=pos_sel_decay)

mixnmatch_lam <- dredge(glb1,rank = "AIC",m.lim = c(0,6)) #6 seems like the maximum terms we can fit safely
av <- model.avg(mixnmatch_lam)
head(mixnmatch_lam, 76)

#Global model call: lm(formula = lambda ~ s41 + s93 + s102 + s160 + s177 + s189 +
#    s261 + s291 + s389 + s477, data = pos_sel_lam)
#---
#Model selection table
#    (Intrc) s102 s160 s177 s189 s261 s291 s389 s41 s477 s93 df logLik  AIC delta weight
#59    8.953         +         +    +    +                    9 -1.117 20.2  0.00  0.017
#187   8.953         +         +    +    +        +           9 -1.117 20.2  0.00  0.017
#315   8.953         +         +    +    +             +      9 -1.117 20.2  0.00  0.017
#443   8.953         +         +    +    +        +    +      9 -1.117 20.2  0.00  0.017
#155   1.458         +         +    +             +           9 -1.117 20.2  0.00  0.017
#281   5.720                   +    +                  +      9 -1.117 20.2  0.00  0.017
#411   1.458         +         +    +             +    +      9 -1.117 20.2  0.00  0.017
#57   13.210                   +    +    +                    9 -1.117 20.2  0.00  0.017
#185  13.210                   +    +    +        +           9 -1.117 20.2  0.00  0.017
#313  13.210                   +    +    +             +      9 -1.117 20.2  0.00  0.017
#441  13.210                   +    +    +        +    +      9 -1.117 20.2  0.00  0.017
#283   8.953         +         +    +                  +      9 -1.117 20.2  0.00  0.017
#153   5.720                   +    +             +           9 -1.117 20.2  0.00  0.017
#409   5.720                   +    +             +    +      9 -1.117 20.2  0.00  0.017
#650   5.883    +              +                  +        + 10 -0.730 21.5  1.22  0.009
#652
#
#use this function to look at each model
summary(eval(getCall(mixnmatch_lam,'59')))


#looking at the two sites are that most commonly present
decay_anova <- lm(lambda ~ s41 + s160 + s189 + s261 + s291 + s477, data = pos_sel_decay) # > 13 appearances
anova(decay_anova)->table5


displaysites(c(41, 189, 261)) -> positivelyselectedThatpredictDecay
positivelyselectedThatpredictDecay
