library(ggplot2)

maindir <- "~/Documents/GitHub/Cypridinidae_EmissionSpectra/"
setwd(maindir)
data <- read.csv(file="Raw Data/expression-kinetics/2014_mammalian_assays.csv", header=TRUE)
#Add column for dilution*substrate
luciferin <- data$dilution*data$substrate
data2 <- cbind(data,luciferin)
data <- data2

#Subset data. the Soro-luc constructs did not have proper signal peptides
tsuluc <- grep("VtLa", data$construct); tsujii <- data[tsuluc,]; 
pmorluc <- grep("PhM", data$construct); morini <- data[pmorluc,]
svuluc <- grep("SVUluc", data$construct); svu <- data[svuluc,]
blank <- grep("Blank", data$construct); blankcells <- data[blank,]
hek <- grep("HEK", data$construct); hekcells <- data[hek,]
combine <- rbind(tsujii, morini, hekcells, blankcells, svu)
subset(combine, log(luciferin) > 6.5) -> maxluc

ggplot(data=maxluc, aes(x=construct, y=log10(light), fill=construct)) + scale_fill_manual(values=c("grey","grey","#0092ff","#0092ff","#0092ff")) + geom_boxplot() + geom_jitter() + xlab("Species") + ylab(expression('log'[10]*'( Counts Per Second)'))