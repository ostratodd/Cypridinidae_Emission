alldata			#This table is constructed in 3_Calculate_parameters.R
alldata -> TableS3
	#Now write all data to text file for supplement
write.table(alldata, file="TableS3.txt", sep="\t")
