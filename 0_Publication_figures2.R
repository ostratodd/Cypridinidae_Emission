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

################
#Table 2 - ANOVA results for 5 ‘mutagenesis sites’ using all data from mutagenesis study and available species’ luciferases
allmutsites <- c(38, 45, 75, 79, 87, 126, 167, 170, 178, 191, 197, 223, 258, 276, 280, 372, 375, 403, 404, 405, 406, 407, 479)
mutsites <- c(38, 178, 375, 404, 405) #Cypridina numbering sites with 3 states
cyptoalign$Aligned[match(allmutsites,cyptoalign$Cypridina_noctiluca_BAD08210)] -> alignment_numbers
anova(lm(lmax ~ X38 * X178 * X375 * X404 * X405, data=cn)) -> table2
write.table(table2, file = "Table2.txt", sep="\t")
table2
print("Tranlastion of Alingment sites to Cypridina site numbers")
allmutsites
alignment_numbers
#rbind(allmutsites, alignment_numbers)

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

#***************************Figure 4 Luciferase Tree and interesting sites
####################
#BiocManager::install("ggtree")
require(ggtree);
library(gridExtra); library(grid);
luc_tree <- read.newick(file="LuciferaseTree_dNds/results/phylogenies/combined_codon.treefile")
root(luc_tree, c("Vargula_hilgendorfii_AAA30332", "Cypridina_noctiluca_BAD08210")) -> luc_tree_root
#CURRENTLY must read positively selected site with code below

data.frame(pos_sel[-1]) -> ps_df; ps_df->ps_df_named;
rownames(ps_df_named) <- pos_sel$sp
#need to order rows the same as tree
ps_plot <- ggtree(luc_tree_root) + geom_tiplab(align=TRUE)

mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 1)),
    colhead = list(fg_params=list(cex = 1)),
    rowhead = list(fg_params=list(cex = 1)))

tt3 <- ttheme_default(
  core=list(bg_params = list(fill = blues9[1:2], col=NA),
            fg_params=list(fontface=1)),
  colhead=list(fg_params=list(col="navyblue", fontface=1)),
  rowhead=list(fg_params=list(col="white")));

quartz("Figure 4", 15, 5) 
tableGrob(ps_df_named, theme=mytheme)->ps_table
grid.arrange(ps_plot, ps_table, ncol=2)



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
#Figure XX - Do + selected sites predict color or kinetics?
library(gridExtra); library(grid); require(janitor);
library(MuMIn);

#Read in luciferase alignment and sites under selection
dat <- read.csv("LuciferaseTree_dNds/results/combined_aa.csv",header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))
#Columns are +1 compared to meme due to sp colum
#Next command alters numbers by 1 to account for this and name species column sp
colnames(dat)[2:ncol(dat)] <- paste("s",seq(1,(ncol(dat)-1)),sep="")
colnames(dat)[1] <- "sp"
require(tidyverse)
dat %>% arrange(sp) -> dat #sort by species name 


#Here is aligned amino acid file used for codon alignment. 
#Read results in csv from meme selection analysis, including positively selected sites
meme <- read.csv("LuciferaseTree_dNds/results/hyphy/lucclade.meme.csv",header=TRUE)
#Pull positively selected sites using vector of output table from meme
pos_sel <- dat[,c("sp",paste("s",meme$Codon, sep=""))]
#decay file has translation of different codes between datasets
decay <- read.csv("Raw Data/expression-kinetics/decay_averages_all_for comparison_with_color.csv",header=TRUE)


#***********color
table1 <- read.table(file="Table1.txt", sep="\t", header=TRUE) #Read again if not executed above
translate <- read.csv("Raw Data/expression-kinetics/translate_lucname_decayname.csv",header=TRUE)
tmpmerge <- merge(translate, table1, by='Species')
lucNcolor <- merge(dat, tmpmerge, by='sp')
pos_sel_col <- lucNcolor[,c("sp",paste("s",meme$Codon, sep=""),"Lmax_Mean", "FWHM_Mean")]
remove_constant(pos_sel_col)->pos_sel_col ##Invariant sites because luc varies in spp without color data 

options(na.action = "na.fail")
#doesn't work inside lm function but can cut and paste manually -- no invariant sites for color
paste(colnames(pos_sel_col)[4:ncol(pos_sel_col)-2], collapse=" + "   )

colorlm_fel <- lm(Lmax_Mean ~ s43 + s209, data=pos_sel_col)


colorlm <- lm(Lmax_Mean ~ s41 + s93 + s102 + s142 + s160 + s177 + s189 + s261 + s285 + s291 + s320 + s389 + s477, data=pos_sel_col)
fwhmlm <- lm(FWHM_Mean ~ s41 + s93 + s102 + s142 + s160 + s177 + s189 + s261 + s285 + s291 + s320 + s389 + s477, data=pos_sel_col)
      
mixnmatch_col <- dredge(	colorlm,rank = "AIC",m.lim = c(0,3)) #3 seems like the maximum terms we can fit safely
head(mixnmatch_col, 12)
av <- model.avg(mixnmatch_col)

#Global model call: lm(formula = Lmax_Mean ~ s41 + s93 + s102 + s142 + s160 + s177 + 
#    s189 + s261 + s285 + s291 + s320 + s389 + s477, data = pos_sel_col)
#---
#Model selection table 
#     (Intrc) s102 s142 s160 s177 s189 s261 s285 s291 s320 s389 s41 s477 s93 df  logLik  AIC delta weight
#51     458.8         +              +    +                                   9  10.990 -4.0  0.00  0.167
#113    458.8                        +    +    +                              9  10.990 -4.0  0.00  0.167
#305    458.8                        +    +              +                    9  10.990 -4.0  0.00  0.167
#82     458.8    +                   +         +                              9  10.990 -4.0  0.00  0.167
#274    458.8    +                   +                   +                    9  10.990 -4.0  0.00  0.167
#20     458.8    +    +              +                                        9  10.990 -4.0  0.00  0.167

summary(eval(getCall(mixnmatch_col,51))) #41
summary(eval(getCall(mixnmatch_col, 113))) #

lmax_anova <- lm(Lmax_Mean ~ s102 + s142 + s189 + s261 + s285 + s320,data=pos_sel_col) # rep'd 3 or more x in best models
anova(lmax_anova)
lmax_anova2 <- lm(Lmax_Mean ~ s93 * s102 * s142 * s189 * s261 * s320, data=pos_sel_col) # rep'd 2 or more x in best models
anova(lmax_anova2)
summary(lmax_anova2)

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


#***********kinetics/decay
### decay ANOVA  ###
#Now read decay data and merge -- decay column called lambda
decay <- read.csv("Raw Data/expression-kinetics/decay_averages_all_for comparison_with_color.csv",header=TRUE)
translate <- read.csv("Raw Data/expression-kinetics/translate_lucname_decayname.csv",header=TRUE)
tmpmerge <- merge(decay, translate, by="PubName")
lucNkinetics <- merge(dat, tmpmerge, by="sp")
pos_sel_lam <- lucNkinetics[,c("sp",paste("s",meme$Codon, sep=""),"lambda")]
remove_constant(pos_sel_lam)->pos_sel_lam ##Invariant sites because luc varies in spp without decay data 

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
# 75 models are within 2 AIC

#use this function to look at each model
summary(eval(getCall(mixnmatch_lam,'59'))) 


#looking at the two sites are that most commonly present
decay_anova <- lm(lambda ~ s160 + s189 + s261 + s291, data = pos_sel_lam) # > 13 appearances
anova(decay_anova)

decay_anova2 <- lm(lambda ~ s41 + s93 + s102 + s160 + s189 + s261 + s291 + s389 + s477, data=pos_sel_lam) # > 12 appearances
anova(decay_anova2)

#############################
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
