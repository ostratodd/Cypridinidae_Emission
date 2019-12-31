maindir <- "~/Documents/GitHub/Cypridinidae_Emission/"
setwd(maindir)

#First read in data sets and complete analyses for figures and tables
source("1_Emission_Functions.R")
source("2_Read_Datasets.R")
source("3_Calculate_parameters.R")
source("4_MutagenesisLuciferase.R")

#*************Main Figures
#Figure 1 is luciferase phylogeny and positively selected (meme) sites
source("Figure1.R")

#Figure 2 is in vitro expression experiments
source("Figure2.R")

#Figure 3 plot lambda max values (color)
source("Figure3.R")
Figure3



#*************************************************Main Tables
#Table 1 - Stats table for in vitro expression

#Table 2 - lmax and fwhm for all species studied
source("Table2.R")

#Table 3 - Sites affecting lmax including mutagenesis data
source("Table3.R")
table3

#Table 4 - positively selected sites correlated with lambda max
source("Table4.R")
table4

#Table 5 - positively selected sites correlated with with decay


#***************************************************Supplemental Figures
#Figure S1 - Luciferase phylogeny with non-luc outgroups
source("FigureS1.R")

#Figure S2 - Plot of all FWHM data
source("FigureS2.R")
FigureS2

#Figure S3 - No correlation between lambda max and decay
source("FigureS3.R")
FigureS3

#**********************************Supplemental Tables
#Supplemental Table S1 is previously published emission spectra (constructed manually outside R)

#Supplemental Table S2 -- Statistics for comparisons during in vitro expression
source("Table")

#Supplemental Table S3- all emission parameter data for each individual organism
source("TableS2.R")

#Supplemental Table S4 - collection and sequencing data

#Supplemental Table S5 - multiple sequence alignment









#***************************** Below here are commands for tables that are not yet separated into their own files



#*************************************************Supplemental Tables


############
#Figure XX - Do positively selected sites predict color or kinetics?
library(gridExtra); library(grid); require(janitor);
library(MuMIn);



#***********color

options(na.action = "na.fail")
#doesn't work inside lm function but can cut and paste manually -- no invariant sites for color
paste(colnames(pos_sel_col)[4:ncol(pos_sel_col)-2], collapse=" + "   )

colorlm_fel <- lm(Lmax_Mean ~ s43 + s209, data=pos_sel_col)



#FWHM
mixnmatch_fw <- dredge(fwhmlm,rank = "AIC",m.lim = c(0,3)) #3 seems like the maximum terms we can fit safely
head(mixnmatch_fw, 12)

#Global model call: lm(formula = FWHM_Mean ~ s41 + s93 + s102 + s142 + s160 + s177 +
#    s189 + s261 + s285 + s291 + s320 + s389 + s477, data = pos_sel_col)
#---
#Model selection table
#     (Intrc) s102 s160 s177 s189 s261 s291 s41 s93 df  logLik  AIC delta weight
#1073   82.52                   +    +        +     10   1.955 16.1  0.00  0.167
#177    84.34                   +    +    +         10   1.955 16.1  0.00  0.167
#57     84.34              +    +    +              10   1.955 16.1  0.00  0.167
#26     84.34    +         +    +                   10   1.955 16.1  0.00  0.167
#146    84.34    +              +         +         10   1.955 16.1  0.00  0.167
#1042   82.52    +              +             +     10   1.955 16.1  0.00  0.167

#4114   88.58    +              +                 +  9  -6.271 30.5 14.45  0.000
#6

fwhm_anova <- lm(FWHM_Mean ~ s102 + s189 + s261 ,data=pos_sel_col) # in 3 or more of top models
anova(fwhm_anova)
fwhm_anova2 <- lm(FWHM_Mean ~ s41 + s102 + s177 + s189 + s261 + s291, data=pos_sel_col) # in 2 or more of top models
anova(fwhm_anova2)
#No positively selected sites predict FWHM

#***********kinetics/decay
### decay ANOVA  ###
#Now read decay data and merge -- decay column called lambda

options(na.action = "na.fail")
#doesn't work inside lm function but can cut and paste manually
#after removing invariant sites, which are 160, 303, 338
paste(colnames(pos_sel_lam)[3:ncol(pos_sel_lam)-1], collapse=" + "   )
#IF sites change copy/paste from result of above command
glb1 <- lm(lambda ~ s41 + s93 + s102 + s160 + s177 + s189 + s261 + s291 + s389 + s477, data=pos_sel_lam)

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
decay_anova <- lm(lambda ~ s160 + s189 + s261 + s291, data = pos_sel_lam) # > 13 appearances
anova(decay_anova)

decay_anova2 <- lm(lambda ~ s41 + s93 + s102 + s160 + s189 + s261 + s291 + s389 + s477, data=pos_sel_lam) # > 12 appearances
anova(decay_anova2)

displaysites(c(41, 102, 189)) -> positivelyselectedThatpredictDecay
positivelyselectedThatpredictDecay
