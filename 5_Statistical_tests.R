### statistical tests for the expression assay data ###
#source publication file first 0_Publications.R
library(lme4); library(rcompanion); library(nlme); library(ggplot2)

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
##############


##########****************************
##HEK mammalian expression statistics###
#this is the same as the data from 0_Publication_figures2.R
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

#t.tests to compare the activity of luciferases expressed in mammalian cells to the HEK cells
mam.pmo <- maxluc[maxluc$construct == "PhM" | maxluc$construct == "HEK",]
t.test(log10(mam.pmo$light)~mam.pmo$construct)
#p-value = 0.001253

mam.svu <- maxluc[maxluc$construct == "SVUluc" | maxluc$construct == "HEK",]
t.test(log10(mam.svu$light)~mam.svu$construct)
#p-value = 4.51e-09

mam.vts <- maxluc[maxluc$construct == "VtLa" | maxluc$construct == "HEK",]
t.test(log10(mam.vts$light)~mam.vts$construct)
#p-value = 4.56e-09

p.adjust(c(0.001253,4.51e-09,4.56e-09),method = "bonf")
#corrected: 3.759e-03, 1.353e-08, 1.368e-08

##########****************************
##Pichia expression statistics###

##mixed model to account for repeated measures amongst samples##
light_data <- read.csv("Raw Data/expression-kinetics/2019_pichia_assays.csv",header=TRUE)
light_data2 <- light_data %>% group_by(Species,type,sample,count_type,Concentration,date) %>% summarise(mean_cps = mean(cps_norm))
light_data_drop <- subset.data.frame(light_data,!is.na(light_data$cps_norm))
light_data_drop$Species <- factor(light_data_drop$Species,levels = c("Pichia","CNO","KHC","SVU","VTS"))
m.rand <- lme(log10(cps_norm) ~ Species*count_type, random = ~1 | sample/Replicate, data = light_data_drop)
m.rand_null <- gls(log10(cps_norm) ~ Species*count_type, data = light_data_drop)
summary(m.rand)
#these are all in comparison to Pichia, after + luciferin as the base intercept. Thus p-values indicate whether the effect is different from that standard.
#thus it seems like generally, CNO (both before and after) is different from Pichia/after
#all species have a significant interaction term, suggesting that each luciferase activity is different from Pichia/after and that changes across before/after measures

anova(m.rand,m.rand_null)
nagelkerke(m.rand,m.rand_null)
plotNormalHistogram(residuals(m.rand))
plot(fitted(m.rand), residuals(m.rand))

##Supplemental figure of the effect sizes (only works with lme4 models)
m.rand_lmer <- lmer(log10(cps_norm) ~ Species*count_type + (1|sample/Replicate),data = light_data_drop)
ggCaterpillar(ranef(m.rand_lmer)) #significant, varying random effects across samples indicate that we should account for it

##t-tests between before/after in each species to test change in light##
cno <- subset.data.frame(light_data2,light_data2$Species == "CNO")
cno <- subset.data.frame(cno,!is.na(cno$mean_cps))
t.test(log10(cno$mean_cps)~cno$count_type)
#p-value = 0.0002923

khc <- subset.data.frame(light_data2,light_data2$Species == "KHC")
t.test(log10(khc$mean_cps)~khc$count_type)
#p-value = 0.02686

svu2 <- subset.data.frame(light_data2,light_data2$Species == "SVU")
t.test(log10(svu2$mean_cps)~svu2$count_type)
#t = 2.6185, df = 11.504, p-value = 0.02313

vts <- subset.data.frame(light_data2,light_data2$Species == "VTS")
t.test(log10(vts$mean_cps)~vts$count_type)
#p-value = 0.005704

p.adjust(c(0.0002923,0.02686,0.02313,0.005704),method = "bonf")
#corrected:0.0011692, 0.1074400, 0.0925200, 0.0228160

