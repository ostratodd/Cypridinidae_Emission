maindir <- "~/Documents/GitHub/Cypridinidae_EmissionSpectra/"
setwd(maindir)

#Test if Vargula tsujii VtL is a functional luciferase
library(ggplot2)
data <- read.csv(file="Raw Data/expression-kinetics/2014_mammalian_assays.csv", header=TRUE)
#Add column for dilution*substrate
luciferin <- data$dilution*data$substrate
data2 <- cbind(data,luciferin)
data <- data2

#HEK -- negative control comparison is untransfected HEK cells
rows.to.keep <- grep("hek", data$genefamily)
hekdata <- data[rows.to.keep,]
ggplot(hekdata, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(date)),size=3) 
#Only include 2013 data 
hek13 <- subset(hekdata, year == 13)
#Only include 2014 data
hek14 <- subset(hekdata, year == 14)
ggplot(hek13, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(date)),size=3) 

#blank -- negative control comparison is untransfected HEK cells
rows.to.keep <- grep("blank", data$genefamily)
blankdata <- data[rows.to.keep,]
ggplot(blankdata, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(date)),size=3) 
#Only include 2013 data 
blank13 <- subset(blankdata, year == 13)
#Only include 2014 data
blank14 <- subset(blankdata, year == 14)
ggplot(blank13, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(date)),size=3) 


#tsujii -- Plot All V. tsujii data
rows.to.keep <- grep("Vargula_tsujii", data$species)
Vtsujiidata <- data[rows.to.keep,]
ggplot(Vtsujiidata, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(construct)),size=4)
#Note that VtLb seems to not be functional. Hypothesis is that clone is problematic.


#Separate 2013 data from 2014 data. There were differences in experimental techniques, including different cell culture media.

#Only include 2013 data 
Vtsujiidata13 <- subset(Vtsujiidata, year == 13)
#Only include 2014 data
Vtsujiidata14 <- subset(Vtsujiidata, year == 14)

ggplot(Vtsujiidata13, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(construct)),size=4)


###Comparison of V tsujii to HEK negative control RAW plot
VtLtest <- rbind(hek13,Vtsujiidata13)
ggplot(VtLtest, aes(x=log(luciferin),y=log(light))) + geom_point(aes(colour=factor(construct)),size=4)

##VtL to HEK aov
mod1 <- aov(light~luciferin*genefamily, data=VtLtest)
summary(mod1)
mod2 <- aov(light~luciferin+genefamily, data=VtLtest)
summary(mod2)
anova(mod1,mod2)

luc <- subset(VtLtest, genefamily == "luciferase")
hek <- subset(VtLtest, genefamily == "hek")

reg1 <- lm(log(light)~log(luciferin), data=luc); summary(reg1)
reg2 <- lm(log(light)~log(luciferin), data=hek); summary(reg2)

plot(log(light)~log(luciferin), data=VtLtest, type='n', xlab="Relative Luciferin Concentration (log)", ylab="Light Production (log CPS)")
points(log(luc$luciferin),log(luc$light), pch=20)
points(log(hek$luciferin),log(hek$light), pch=1)
abline(reg1, lty=1)
abline(reg2, lty=2)
legend("bottomright", c("VtLa","HEK293"),lty=c(1,2),pch=c(20,1))

#Make Boxplot of 2014 data
#Remove VtLb data
VtLa14 <- subset(Vtsujiidata14, construct == "VtLa")
VtLbox <- rbind(VtLa14,hek14,blank14)
VtLbox <- droplevels(VtLbox)
qplot(genefamily, log(light), data=VtLbox, geom="boxplot")
boxplot(log(light)~genefamily, data=VtLbox, ylab="Light Production (log CPS)")

#Plot both main plots together
#par(mfrow=c(1,2))
layout(matrix(c(1,2),1,byrow=TRUE),widths=c(1,2))
boxplot(log(light)~genefamily, data=VtLbox, xlab="Fixed Luciferin Conc.", ylab="Light Production (log CPS)")

plot(log(light)~log(luciferin), data=VtLtest, type='n', xlab="Luciferin Conc.", ylab="")
points(log(luc$luciferin),log(luc$light), pch=20)
points(log(hek$luciferin),log(hek$light), pch=1)
abline(reg1, lty=1)
abline(reg2, lty=2)
legend("bottomright", c("VtLa","HEK293"),lty=c(1,2),pch=c(20,1))
