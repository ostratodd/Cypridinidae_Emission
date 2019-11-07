#Mac Laptop
setwd("~/Documents/GitHub/Cypridinidae_EmissionSpectra/Luciferase/")
#test with previously published analysis by Yokoyama
read.csv(file="YokoyamaOpsin.csv") -> opsin

#omnibus test two way
anova(lm(lmax ~ S180 * S197 * S277* S285 * S308, data=opsin))

#oneway
(lm(lmax ~ pigmenttype, data=opsin))

###Cypridina mutations and lamda-max

read.csv(file="Cypridina only.csv") -> cn
#Calculate number of factors (different amino acids) per site
sapply(cn, function(x) length(unique(x))) -> aapersite

#Choose sites (there are 5) with >2 different factors (amino acids at that site)
anova(lm(lmax ~ X38 * X178 * X375 * X404 * X405, data=cn))
plot.design(lmax ~  X38 * X178 * X375 * X404 * X405, data=cn)

#now recalculate adding non-cypridina species
read.csv(file="Cypridina mutagenesis.csv") -> cnplus
anova(lm(lmax ~ X38 * X178 * X375 * X404 * X405, data=cnplus))
plot.design(lmax ~  X38 * X178 * X375 * X404 * X405, data=cnplus)

plot(cnplus$Nsubs, cnplus$lmax)

#Now calculate with non-cypridina only -- remove 405 which is constant
read.csv(file="NonCypridina only.csv") -> nc
anova(lm(lmax ~ X38 * X178 * X191, data=nc))
plot.design(lmax ~  X38 * X178 * X375 * X404, data=nc)

#There of the 23 sites mutated in patent application, there are 4 sites conserved in our 3 Photeros but also different in non-Photeros

read.csv(file="Cypridina mutagenesis.csv") -> cnplus
anova(lm(lmax ~ X178 * X191 * X223 * X403, data=cnplus))
plot.design(lmax ~  X178 * X191 * X223 * X403, data=cnplus)