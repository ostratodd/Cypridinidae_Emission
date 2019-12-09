maindir <- "~/Documents/GitHub/Cypridinidae_EmissionSpectra/"
setwd(maindir)

source("1_Emission_Functions.R")
source("2_Read_Datasets.R")
source("3_Calculate_parameters.R")
source("4_MutagenesisLuciferase.R")

#*************************************************Main Tables
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
	#now write to text file
write.table(table1, file = "Table1.txt", sep="\t")
#graphical table
library(gridExtra); library(grid);
tt3 <- ttheme_default(
  core=list(bg_params = list(fill = blues9[1:2], col=NA),
            fg_params=list(fontface=1)),
  colhead=list(fg_params=list(col="navyblue", fontface=1)),
  rowhead=list(fg_params=list(col="white")));
grid.arrange(tableGrob(table1, theme=tt3))

#Table 2 - ANOVA results for 5 ‘mutagenesis sites’ using all data from mutagenesis study and available species’ luciferases
anova(lm(lmax ~ X38 * X178 * X375 * X404 * X405, data=cn)) -> table2
write.table(table2, file = "Table2.txt", sep="\t")

#*************************************************Main Figures
################
#Figure 1 - Box plot of all lmax of species
require(ggplot2)
clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
quartz("Figure 1", 13, 3)
fig1 <- qplot(abbreviation, sgMax, data=clean, geom=c("boxplot", "jitter")) + ylab("Lambda max (nm)") + xlab("Species")
	#highlight data with rectangles
	#add green blue shading w/ rectangle (geom_rect). This uses y-co-ordinates as min and max of photeros (green) and non-photeros (blue). facet_grid sorts into panels by collection locality
t2.rect1 <- data.frame (xmin=-Inf, xmax=Inf, ymin=463.97, ymax=471.14)
t1.rect1 <- data.frame (xmin=-Inf, xmax=Inf, ymin=454.27, ymax=463.42)

fig1 + geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) + geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.07, inherit.aes = FALSE)  + facet_grid(cols = vars(locality), scales = "free_x", switch = "x", space = "free")

##############
#Figure 2 - Box plot of all fwhm of species
require(ggplot2)
clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")

quartz("Figure 2", 13, 3)
fig2 <- qplot(abbreviation, sgfwhm, data=clean, geom=c("boxplot", "jitter")) + ylab("FWHM (nm)") + xlab("Species")
fig2 + facet_grid(cols = vars(locality), scales = "free_x", switch = "x", space = "free")

##############
#Figure 3 - Expression of luciferase candidates in Pichia yeast
library(ggplot2); library(reshape2); library(dplyr); library(tidyverse); library(gridExtra)
setwd(maindir)
light_data <- read.csv("Raw Data/expression-kinetics/2019_pichia_assays.csv",header=TRUE)
light_data2 <- light_data %>% group_by(Species,type,sample,count_type,Concentration,date) %>% summarise(mean_cps = mean(cps_norm))

#Pichia yeast expression plot###
sp <- c("Pichia","CNO","KHC","SVU","VTS")
light_data2$Species <- factor(light_data2$Species,levels=sp)
Figure3b <- ggplot(data=light_data2, aes(x = Species, y = log10(mean_cps),fill=count_type)) + geom_boxplot() +
  geom_point(position = position_jitterdodge()) + scale_discrete_manual(aes(x=sp)) +
  xlab("Species") + ylab(expression('log'[10]*'( CPS / Total Protein Conc. )')) +
  scale_fill_manual(values = c("#0092ff","grey"),name="Measures",labels=c("After + luciferin","Before")) +
  scale_x_discrete(labels=c("Pichia","C_noc","K_has","M_SVU","V_tsu")) +
  theme(legend.position = "none")

	##############
	#Mammalian expression
mam_data <- read.csv(file="Raw Data/expression-kinetics/2014_mammalian_assays.csv", header=TRUE)
#Add column for dilution*substrate
luciferin <- mam_data$dilution*mam_data$substrate
mam_data2 <- cbind(mam_data,luciferin)
mam_data <- mam_data2

#Subset data. the Soro-luc constructs did not have proper signal peptides
	tsuluc <- grep("VtLa", mam_data$construct); tsujii <- mam_data[tsuluc,];
	pmorluc <- grep("PhM", mam_data$construct); morini <- mam_data[pmorluc,]
	svuluc <- grep("SVUluc", mam_data$construct); svu <- mam_data[svuluc,]
	blank <- grep("Blank", mam_data$construct); blankcells <- mam_data[blank,]
	hek <- grep("HEK", mam_data$construct); hekcells <- mam_data[hek,]
	mam_combine <- rbind(tsujii, morini, hekcells, blankcells, svu)
	subset(mam_combine, log(luciferin) > 6.5) -> maxluc

Figure3a <- ggplot(data=maxluc, aes(x=construct, y=log10(light), fill=construct)) + scale_fill_manual(values=c("grey","grey","#0092ff","#0092ff","#0092ff")) + geom_boxplot() + geom_jitter() + xlab("Species") + ylab(expression('log'[10]*'( Counts Per Second)')) +
  scale_x_discrete(labels = c("Blank", "HEK", "P_mor", "M_SVU", "V_tsu")) +
  theme(legend.position = "none")

grid.arrange(Figure3a, Figure3b, nrow = 1) -> Figure3


#*************************************************Supplemental Tables
#Supplemental Table S1 is previously published emission spectra
#Supplemental Table S2 is collecting information
#Supplemental Table S3 -- all emission parameter data for each individual organism
alldata
	#Now write all data to text file for supplement
write.table(alldata, file="TableS3.txt", sep="\t")
#Supplemental Table S4
anova(lm(lmax ~ X38 * X178 * X191, data=nc)) -> tableS4
write.table(tableS4, file = "TableS4.txt", sep="\t")
#Supplemental Table S5 - Full table of mutant and natural luciferases and corresponding emission spectra
cnplus->tableS5
write.table(tableS5, file = "TableS5.txt", sep="\t")

#***************************************************Supplemental Figures
#Supplemental Figure S1
require(ape)
require(phytools)
read.tree(file="./LuciferaseTree_dNds/results/phylogenies/all_aa.treefile") -> lucplus_tree
midpoint.root(lucplus_tree) -> lucplus_tree_r
ladderize(lucplus_tree_r)-> lprl
plot(lprl, show.tip.label=TRUE, cex=.8, x.lim=2)

#Supplemental Figure S2
#Created with prettyplot

############
#Supplemental Figure S3
### decay ANOVA for luciferase function paper ###
# using codon-aligned translations from AliView given THO's dNds output

dat <- read.csv("Raw Data/expression-kinetics/clipboard-alignment_4671897282832913936.translated.csv",header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))
#sites from positive selection scan are:
# 93
# 115
# 142
# 152
# 160
# 189
# 261
# 285
# 320
# 371
# 389
# 477
# 506
# 581

#sites from color data converted to this alignment are, converted to CNO alignment:
# 38 --> 64
# 87 --> 115
# 178 --> 207
# 375 --> 406
# 404 --> 435
# 405 --> 436 ##invariant in decay dataset

#naming all columns properly
colnames(dat)[1] <- "sp"
colnames(dat)[589] <- "decay"
<<<<<<< HEAD

#sites that are invariant in the list above are:
# 64
# 115
# 142
# 285
# 320
# 371
# 406
# 506
library(gridExtra); library(grid);
dat[,154]

dat2 <- subset.data.frame(dat,dat$sp != "Vargula_tsujii_sequenced")

#add one for each to account for species names in first column
#Positive selection sites -- 390 causes error
m.1 <- lm(decay~ (V94 + V116 + V143 + V153 + V161 + V190 + V262 + V286 + V321 + V372 + V478 + V507 + V582),data=dat2)
anova(m.1)

m.2 <- lm(decay ~ V94 + V153 + V190 + V478,data=dat2)
anova(m.2)
=======
colnames(dat)[2:588] <- paste("V",seq(1,587),sep="")
dat2 <- dat[-8,]
dat2 <- dat2[,c("sp","V93","V115","V142","V152","V160","V189","V261","V285",
               "V320","V371","V389","V477","V506","V581","decay")]

library(MuMIn)
options(na.action = "na.fail")
glb1 <- lm(decay ~ V93 + V142 + V152 + V160 + V189 + V261 + V285 + V320 + V371 + V389 + V477 + V506 + V581, data = dat2)
mixnmatch <- dredge(glb1,rank = "AIC",m.lim = c(0,6)) #6 seems like the maximum terms we can fit safely
av <- model.avg(mixnmatch)

#model averaging approach produces these top 6 equivalent models:
#       (Intrc) V142 V152 V160 V189 V261 V285 V320 V371 V389 V477 V506 V581 V93 df logLik  AIC  delta weight
# 796   8.8620    +    +         +    +                   +    +               12  15.512  -7.0  0.00  0.166
# 827   8.8620         +         +    +    +              +    +               12  15.512  -7.0  0.00  0.166
# 859   8.8620         +         +    +         +         +    +               12  15.512  -7.0  0.00  0.166
# 2332 15.7400    +    +         +    +                   +              +     12  15.512  -7.0  0.00  0.166
# 2363 15.7400         +         +    +    +              +              +     12  15.512  -7.0  0.00  0.166
# 2395 15.7400         +         +    +         +         +              +     12  15.512  -7.0  0.00  0.166

#use this function to look at each model
summary(eval(getCall(mixnmatch,796))) #189, 581
summary(eval(getCall(mixnmatch,827))) #371, 477, 581
summary(eval(getCall(mixnmatch,859))) #371, 581
summary(eval(getCall(mixnmatch,2332))) #506
summary(eval(getCall(mixnmatch,2363))) #371
summary(eval(getCall(mixnmatch,2395))) #506

#looking at the four sites are that ALWAYS present
m_always <- lm(decay ~ V152 + V189 + V261 + V389,data=dat2)
anova(m_always)

m_some <- lm(decay ~ V142 + V285 + V320 + V477 + V581,data=dat2)
anova(m_some)

m_never <- lm(decay ~ V93 + V160 + V371 + V506,data=dat2)
anova(m_never)

#looking at a model that has all sites that appeared significant in at least one of the top models / comparisons above
mixer <- lm(decay ~ V189 + V371 + V477 + V506 + V581,data=dat2)
#site identity is correlated b/c limited sequence diversity
>>>>>>> 22d14059627d67655e952a4e3f935ac2d3d20bee

## plot of lamda max and decay
table1 <- read.table(file="Table1.txt", sep="\t", header=TRUE) #Read again if not executed above
decay <- read.csv("Raw Data/expression-kinetics/decay_averages_all_for comparison_with_color.csv",header=TRUE)
col_dec <- merge(decay,table1,by="Species")
head(col_dec)
plot(Lmax_Mean~lambda,data=col_dec)
cor.test(col_dec$lambda,col_dec$Lmax_Mean)
cor.test(col_dec$lambda,col_dec$FWHM_Mean)
library(ggplot2)

figS <- ggplot(data=col_dec,aes(x=lambda,y=Lmax_Mean)) + geom_point(aes(color=genus,shape=country),size=3) +
  xlab("Ave. decay constant per species") + ylab("Ave. peak emission per species")


###******GGplot custom function to plot lme4 catepillar plot (confidence intervals of random effects)
##credit to davebraze for this function: https://rdrr.io/github/davebraze/FDB1/man/ggCaterpillar.html
##stack exchange resources: https://stackoverflow.com/questions/13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-how-to-mak

ggCaterpillar <- function(re, QQ=FALSE, likeDotplot=TRUE, detailedFacetLabs = TRUE) {
  f <- function(x, nm = "ranef plot") {
    pv   <- attr(x, "postVar")
    cols <- 1:(dim(pv)[1])
    se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
    ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
    pDf  <- data.frame(y=unlist(x)[ord],
                       ci=1.96*se[ord],
                       nQQ=rep(stats::qnorm(stats::ppoints(nrow(x))), ncol(x)),
                       ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                       ind=gl(ncol(x), nrow(x), labels=names(x)))

    if(detailedFacetLabs){
      pDf$ind <- ifelse(grepl("(Intercept)", pDf$ind), "intercept adjustment", paste0("slope adj: ", pDf$ind))
    }

    if(QQ) {  ## normal QQ-plot
      p <- ggplot(pDf, aes_string(x="nQQ", y="y"))
      p <- p + facet_wrap(~ ind, scales="free")
      p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
    } else {  ## caterpillar dotplot
      p <- ggplot(pDf, aes_string(x="ID", y="y")) + coord_flip()
      if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
        p <- p + facet_wrap(~ ind)
      } else {           ## different scales for random effects
        p <- p + facet_grid(ind ~ ., scales="free_y")
      }
      p <- p + xlab(nm) + ylab("Random effects")
      scale <- 12-log(length(levels(pDf$ID)),2)
      p <- p + theme(axis.text.y = element_text(size=scale))
    }

    p <- p + theme(legend.position="none")
    # p <- p + labs(title= nm)
    p <- p + geom_hline(yintercept=0, lwd = I(7/12), colour = I(grDevices::hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
    p <- p + geom_errorbar(aes_string(ymin="y - ci", ymax="y + ci"), width=0, colour="black")
    p <- p + geom_point(aes())
    return(p)
  }

  #   lapply(re, f) # original
  lapply(seq_along(re), function(y, n, i) { f(y[[i]], n[[i]]) }, y=re, n=names(re)) # adds plot names
}
