#Figure 2 - In vitro expression of exemplars supports luciferase function
library(ggplot2); library(reshape2); library(dplyr); library(tidyverse); library(gridExtra)

light_data <- read.csv("Raw Data/expression-kinetics/2019_pichia_assays.csv",header=TRUE)
light_data2 <- light_data %>% group_by(Species,type,sample,count_type,Concentration,date) %>% summarise(mean_cps = mean(cps_norm))

#Pichia yeast expression plot###
sp <- c("Pichia","CNO","KHC","SVU","VTS")
light_data2$Species <- factor(light_data2$Species,levels=sp)
Figure2b <- ggplot(data=light_data2, aes(x = Species, y = log10(mean_cps),fill=count_type)) + geom_boxplot() +
  geom_point(position = position_jitterdodge()) + scale_discrete_manual(aes(x=sp)) +
  xlab("Species") + ylab(expression('log'[10]*'(CPS / Total Protein Conc. )')) +
  scale_fill_manual(values = c("#0092ff","white"),name="Measures",labels=c("After + luciferin","Before")) +
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

quartz("Figure 2", 9, 3)

Figure2a <- ggplot(data=maxluc, aes(x=construct, y=log10(light), fill=construct)) + scale_fill_manual(values=c("white","white","#0092ff","#0092ff","#0092ff")) + geom_boxplot() + geom_jitter() + xlab("Species") + ylab(expression('log'[10]*' (Counts Per Second)')) +
  scale_x_discrete(labels = c("Blank", "HEK", "P_mor", "M_SVU", "V_tsu")) +
  theme(legend.position = "none")

grid.arrange(Figure2a, Figure2b, nrow = 1)
