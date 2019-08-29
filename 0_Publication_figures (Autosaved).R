maindir <- "~/Documents/GitHub/Cypridinidae_EmissionSpectra/"
setwd(maindir)

source("1_Emission_Functions.R")
source("2_Read_Datasets.R")
source("3_Calculate_parameters.R")


#Table 1 - lmax and fwhm for all species studied

#First remove runs where signal:noise is too low
clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
#calculate means, SDs, and counts, grouped by abbreviation
aggregate(clean[, 9:10], list(clean$abbreviation), FUN=mean) -> means
aggregate(clean[, 9:10], list(clean $abbreviation), FUN=sd) -> sds
aggregate(clean[, 9:10], list(clean $abbreviation), FUN=length) -> counts
#Now construct data frame as table for above information
cbind(subset(means, select="Group.1"), counts$sgMax, signif(means$sgMax,4), signif(sds$sgMax,2), signif(means$sgfwhm,4), signif(sds$sgfwhm,2)) -> table1
colnames(table1)[c(1,2,3,4,5,6)] <- c("Species", "N", "Lmax_Mean", "Lmax_SD", "FWHM_Mean", "FWHM_SD")
table1

#Supplemental Table -- all data
alldata

#Figure 1 - Box plot of all lmax of species
require(ggplot2)
clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
#quartz("Figure 1", 13, 3)
fig1 <- qplot(abbreviation, sgMax, data=clean, geom=c("boxplot", "jitter"))
#highlight data with rectanbles
t2.rect1 <- data.frame (xmin=-Inf, xmax=Inf, ymin=463.97, ymax=471.14)
t1.rect1 <- data.frame (xmin=-Inf, xmax=Inf, ymin=454.27, ymax=463.42)

fig1 + geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) + geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.07, inherit.aes = FALSE) + facet_wrap(~locality, scales="free_x")


#Figure 2 - Box plot of all fwhm of species
require(ggplot2)
#a <- subset(alldata, error < 0.02 | replicate=="Cn" | replicate=="Pg")
clean <- subset(alldata, error < 0.02)

quartz("Figure 2", 13, 3)
p <- qplot(abbreviation, sgfwhm, data=clean, geom=c("boxplot", "jitter"))
p 

#*********************************************************************************
#Experimental plots
#Show relationship between lmax and fwhm is different between Photeros and others
par(mfrow=c(1,2))
plot(means$sgMax, means$sgfwhm, xlim=c(460, 462), title(main="Non-Photeros max vs width"))
plot(means$sgMax, means$sgfwhm, xlim=c(465, 469), title(main="Photeros max vs width"))

#Check for correlation between lmax and fwhm************************************************************************************
#Check for correlation of fwhm and lmax in Photeros
subset(clean, genus =="Photeros") -> photeros
aggregate(photeros[, 9:10], list(photeros$abbreviation), FUN=mean) -> photerosmeans
plot(photerosmeans$sgMax, photerosmeans$sgfwhm, xlim=c(465, 472), title(main="Photeros max vs width"))
lmPhoteros = lm(sgMax~sgfwhm, data = photerosmeans)
summary(lmPhoteros)

	#Published data point - so using very different methods -- is a strong outlier
subset(clean, genus =="Photeros" & locality != "Published") -> photeros
aggregate(photeros[, 9:10], list(photeros$abbreviation), FUN=mean) -> photerosmeans
plot(photerosmeans$sgMax, photerosmeans$sgfwhm, xlim=c(465, 470), title(main="Photeros max vs width"))
lmPhoteros = lm(sgMax~sgfwhm, data = photerosmeans)
summary(lmPhoteros)

#Check for correlation of fwhm and lmax in non-Photeros
subset(clean, genus !="Photeros") -> nonphoteros
aggregate(nonphoteros[, 9:10], list(nonphoteros$abbreviation), FUN=mean) -> nonphoterosmeans
plot(nonphoterosmeans$sgMax, nonphoterosmeans$sgfwhm, xlim=c(450, 462), title(main="Non-Photeros max vs width"))
lmNonPhoteros = lm(sgMax~sgfwhm, data = nonphoterosmeans)
summary(lmNonPhoteros)

	#Again, published data point is outlier
subset(clean, genus !="Photeros" & locality != "Published") -> nonphoteros
aggregate(nonphoteros[, 9:10], list(nonphoteros$abbreviation), FUN=mean) -> nonphoterosmeans
plot(nonphoterosmeans$sgMax, nonphoterosmeans$sgfwhm, xlim=c(458, 462), title(main="Non-Photeros max vs width"))
lmNonPhoteros = lm(sgMax~sgfwhm, data = nonphoterosmeans)
summary(lmNonPhoteros)

