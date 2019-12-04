maindir <- "~/Documents/GitHub/Cypridinidae_EmissionSpectra/"
setwd(maindir)

#Test if Vargula tsujii VtL is a functional luciferase
library(ggplot2)
data <- read.csv(file="Raw Data/expression-kinetics/2014_mammalian_assays.csv", header=TRUE)
#Add column for dilution*substrate
luciferin <- data$dilution*data$substrate
data2 <- cbind(data,luciferin)
data <- data2

tsuluc <- grep("VtLa", data$construct)
tsujii <- data[tsuluc,]
ggplot(tsujii, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(year)),size=4)

pmorluc <- grep("PhM", data$construct)
morini <- data[pmorluc,]

svuluc <- grep("SVUluc", data$construct)
svu <- data[svuluc,]


blank <- grep("Blank", data$construct)
blankcells <- data[blank,]

hek <- grep("HEK", data$construct)
hekcells <- data[hek,]



combine <- rbind(tsujii, morini, hekcells, blankcells, svu)
ggplot(combine, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(construct)),size=4)

subset(combine, log(luciferin) > 6.5) -> maxluc
qplot(construct, log(light), data=maxluc, geom="boxplot")
boxplot(data=maxluc, log(light)~construct,  ylab="Light Production (log CPS)")

