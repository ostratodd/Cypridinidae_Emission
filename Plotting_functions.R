#A function to plot raw data (datafile), selecting just one sampling point (column). Note column
#1 is wavelength.
#This may be used to check for saturation of sensor visually
plot1sample <- function(datafile, column) {
	rawdata <- read.table(datafile);
	plot(rawdata[,1], rawdata[,column])
}
	
	
	
	
	
	
