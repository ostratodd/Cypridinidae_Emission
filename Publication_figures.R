maindir <- "~/Desktop/Color_Variation/"
setwd(maindir)


source("1_Emission_Functions.R")
setwd(maindir)
source("2_Read_Datasets.R")
setwd(maindir)
source("3_Calculate_parameters.R")


#Table 1 - lmax and fwhm for all species studied

#First remove runs where signal:noise is too low
#a <- subset(alldata, x < 0.02 | replicate=="Cn" | replicate=="Pg")
a <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
aggregate(a[, 7:8], list(a$abbreviation), mean) -> means
means

#Supplemental Table -- all data
alldata

#Figure X - Box plot of all lmax of species
require(ggplot2)
#a <- subset(alldata, error < 0.02 | replicate=="Cn" | replicate=="Pg")
a <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")

p <- qplot(abbreviation, sgMax, data=a, geom=c("boxplot", "jitter"))
p


#Figure X - Box plot of all fwhm of species
require(ggplot2)
#a <- subset(alldata, error < 0.02 | replicate=="Cn" | replicate=="Pg")
a <- subset(alldata, error < 0.02)

p <- qplot(abbreviation, sgfwhm, data=a, geom=c("boxplot", "jitter"))
p 