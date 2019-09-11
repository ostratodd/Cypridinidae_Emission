#Mac Laptop
setwd("~/Documents/GitHub/Cypridinidae_EmissionSpectra/Luciferase/")
#test with previously published analysis by Yokoyama
read.csv(file="YokoyamaOpsin.csv") -> opsin

#omnibus test two way
anova(lm(lmax ~ S180 * S197 * S277* S285 * S308, data=opsin))

#oneway
(lm(lmax ~ pigmenttype, data=opsin))

###Cypridina mutations and lamda-max

read.csv(file="Cypridina mutagenesis.csv") -> cn

anova(lm(lmax ~ X38 * X75 * X178 * X223 * X280 * X375 * X403, data=cn))

plot.design(lmax ~  X38 * X75 * X178 * X223 * X280 * X375 * X403, data=cn)


plot(cn$Nsubs, cn$lmax)
average()