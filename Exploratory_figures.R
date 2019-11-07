#Must execute 0_Publication before using this
#These are plots used to understand the data that did not make the final publication

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

#*************************************************************
#Does error correlate with distance from mean of that species?
newdata <- subset(alldata, locality != "Published")

#first calculate means using newdata 
aggregate(newdata[, 9:10], list(newdata$abbreviation), FUN=mean) -> means
aggregate(newdata[, 9:10], list(newdata $abbreviation), FUN=sd) -> sds
aggregate(newdata[, 9:10], list(newdata $abbreviation), FUN=length) -> counts
#Now construct data frame as table for above information
cbind(subset(means, select="Group.1"), counts$sgMax, signif(means$sgMax,4), signif(sds$sgMax,2), signif(means$sgfwhm,4), signif(sds$sgfwhm,2)) -> tableS3
colnames(tableS3)[c(1,2,3,4,5,6)] <- c("Species", "N", "Lmax_Mean", "Lmax_SD", "FWHM_Mean", "FWHM_SD")
tableS3
#now compare each observation to mean

errortable <- matrix(ncol=4, nrow=nrow(newdata))
#Here calculate lmax
for (row in 1:nrow(newdata)) {
    lmax <- newdata[row, "sgMax"]
    error <- newdata[row, "error"]
    abbr <- toString(newdata[row,"abbreviation"])
    #find mean of current species value in Table S3 calculated above
    tableS3[which(tableS3 == abbr), 3]
    #create df    
    errortable[row,1] <- lmax
    errortable[row,2] <- error
    tableS3[which(tableS3 == abbr), 3] -> mean
    errortable[row,3] <- mean
    errortable[row,4] <-   abs(lmax-mean)
}
errortable
plot(errortable[,2], errortable[,4])

#Here calculate fwhm
fwerrortable <- matrix(ncol=4, nrow=nrow(newdata))
for (row in 1:nrow(newdata)) {
    fwhm <- newdata[row, "sgfwhm"]
    error <- newdata[row, "error"]
    abbr <- toString(newdata[row,"abbreviation"])
    tableS3[which(tableS3 == abbr), 2] -> nreps
    #find mean of current species value in Table S3 calculated above
    tableS3[which(tableS3 == abbr), 5] #this takes the data from tableS3 which is means and SD by species 5 column is FWHM
    #create df  
    if(nreps > 4){  
    	fwerrortable[row,1] <- fwhm
    	fwerrortable[row,2] <- error
    	tableS3[which(tableS3 == abbr), 5] -> mean
    	fwerrortable[row,3] <- (mean)
    	fwerrortable[row,4] <-  (abs(fwhm-mean))
    }
}
fwerrortable
plot((fwerrortable[,2]),(fwerrortable[,4]))
