maindir <- "~/Documents/GitHub/Cypridinidae_EmissionSpectra/"
setwd(maindir)

source("1_Emission_Functions.R")
source("2_Read_Datasets.R")
source("3_Calculate_parameters.R")


#Table 1 - lmax and fwhm for all species studied

#First remove runs where signal:noise is too low
#a <- subset(alldata, x < 0.02 | replicate=="Cn" | replicate=="Pg")
a <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
#calculate means, SDs, and counts, grouped by abbreviation
aggregate(a[, 7:8], list(a$abbreviation), FUN=mean) -> means
aggregate(a[, 7:8], list(a$abbreviation), FUN=sd) -> sds
aggregate(a[, 7:8], list(a$abbreviation), FUN=length) -> counts
#Now construct data frame as table for above information
cbind(subset(means, select="Group.1"), counts$sgMax, signif(means$sgMax,4), signif(sds$sgMax,2), signif(means$sgfwhm,4), signif(sds$sgfwhm,2)) -> table1
colnames(table1)[c(1,2,3,4,5,6)] <- c("Species", "N", "Lmax_Mean", "Lmax_SD", "FWHM_Mean", "FWHM_SD")
table1

#Supplemental Table -- all data
alldata

#Figure 1 - Box plot of all lmax of species
require(ggplot2)
#a <- subset(alldata, error < 0.02 | replicate=="Cn" | replicate=="Pg")
a <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")

fig1 <- qplot(abbreviation, sgMax, data=a, geom=c("boxplot", "jitter"))
fig1


#Figure X - Box plot of all fwhm of species
require(ggplot2)
#a <- subset(alldata, error < 0.02 | replicate=="Cn" | replicate=="Pg")
a <- subset(alldata, error < 0.02)

p <- qplot(abbreviation, sgfwhm, data=a, geom=c("boxplot", "jitter"))
p 