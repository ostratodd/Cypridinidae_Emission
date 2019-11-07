require(signal)
fullrow <- function(abbreviation, locality, genus, species, replicate, df, sex, preservation, source) {
	#This function creates a full data row adding all metadata
	#and calculating values from spectra, including lmax, fwhm, and error (signal:noise)
	sgsmooth <- 35 #This is the smoothing value for SG smoothing. Doesn't seem to affect results at all
	round(error(df,500),digits=4)->error
	sgfwhm(df,sgsmooth,FALSE)->fwhm
	sgMax(df, sgsmooth)-> lmax
	data.frame(abbreviation, locality, genus, species, replicate, sex, preservation, source, lmax,fwhm,error) -> t1
	#merge(round(error(df,500),digits=4), t4) -> t5
}
graphEmission <- function(datafile, backgroundfile, calibrationfile, cols, plot) {
	rawdata <- read.table(datafile);
	background <-  read.table(backgroundfile);
	calibration <- read.table(calibrationfile);
	
	#USAGE
	#To use this, call the function with datafile (raw), backgroundfile (this is any background, like dim lights etc), calibrationfile 
	#(this is from the blackbody radiation used to standardize), plot (TRUE or FALSE whether to plot or just read data)
	
	
	#***********************Generally not needed, except for early data sets
	##In earlier setups, we learned we Need to extrapolate calibration wavelengths to data wavelengths
	##This is because black body radiation was done on a different day
	##than luminous data collection and the ccd was not aligned precisely the same
	##way each time. This leads to different specific wavelenghts being read
	##Fit high order polynomial to correction factor data. Then interpolate values for data collection wavelengths
	#p20 <- lm( calibration$V2~poly(calibration$V1,20) )
	#xx <- seq(min(calibration$V1),max(calibration$V1), length.out=1024)
	#intercalibration<-data.frame(wavelength=xx, cf=predict(p20, data.frame(x=xx)))
	##Put all data into single data frame
	#fulldata <- cbind(background$V2, intercalibration$cf, rawdata)

	fulldata <- cbind(background$V2, calibration$V2, rawdata)
		

	colnames(fulldata) <- c("background", "calibration", "wavelength");

	#Subtract background from each spectrum
	subtracted <- fulldata[4:(cols+3)]-fulldata$background;
	#Add up all the replicates	
	sumreplicates <- apply(subtracted[,c(1,cols)], 1, sum);
	#multiply by correction factor determined using black body radiator
	calibrated <- sumreplicates * fulldata$calibration

	#Normalize to 1 
    sumreplicates <- sumreplicates/max(sumreplicates, na.rm=T)
    calibrated <- calibrated/max(calibrated, na.rm=T)


	uncalframe <- data.frame(wavelength = fulldata$wavelength, sum=sumreplicates);
	finalframe <- data.frame(wavelength = fulldata$wavelength, sum=calibrated);
	if(plot==TRUE){
		#Uncomment two lines below to plot uncalibrated data
		#plot(uncalframe, xlim=c(400,600), ylim=c(0,1));
		#par(new=TRUE);
		plot(finalframe, xlim=c(400,600), ylim=c(0,1), col="blue");
	}
	#return background-subtracted, summed (over sampling points), and calibrated (using black body radiator) values in data frame by wavelength
	return(finalframe);

}
sgMax <- function(df, smooth) {
#This uses Sovitzky Golay smoothing to smooth sampling curve
#df is a data frame of the corrected data, read in with graphEmission
#smooth is the smoothing value of the SG function
#The function returns the lambda max
	sgdf <- sgolayfilt(df$sum, 3, smooth)
	df2 <- data.frame(wavelength=df$wavelength,sum=sgdf)
	return(data.frame(sgMax=mean(df2[which(df2$sum == max(df2$sum)), ]$wavelength)) )
	
}
sgfwhm <- function(df, smooth, plot) {
#This uses Sovitzky Golay smoothing to smooth sampling curve
#df is a data frame of the corrected data, read in with graphEmission
#smooth is the smoothing value of the SG function
#The function returns the full width half max

	 sgdf <- sgolayfilt(df$sum, 3, smooth)
	 d <- data.frame(wavelength=df$wavelength,sum=sgdf)

	xmax <- d$wavelength[d$sum==max(d$sum)]

	x1 <- d$wavelength[d$wavelength < xmax][which.min(abs(d$sum[d$wavelength < xmax]-max(d$sum)/2))]
	x2 <- d$wavelength[d$wavelength > xmax][which.min(abs(d$sum[d$wavelength > xmax]-max(d$sum)/2))]
	if(plot==TRUE){
		plot(d, type="l") #Can plot with FWHM points
		points(c(x1, x2), c(d$sum[d$wavelength ==x1], d$sum[d$wavelength ==x2]), col="red")
	}

	fwhm <- x2-x1
	data.frame(sgfwhm=fwhm) -> sgfwhm
	return(sgfwhm)
}
error <- function(df, rows) {
	#This is a metric of signal to noise
	#df is data frame with wavelength, sum
	#n is number of highest rows (by sum, which adds all sampled time points together) to average for wavelength
	#the data are ordered, and the "rows" decides how many of the lowest datapoints represent error
	#averaging over them, instead of simply taking the lowest value as the background
	sd(head(df[order(df$sum, decreasing=F),], n=rows)$sum ) -> error
	return(error)
}