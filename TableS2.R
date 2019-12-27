alldata			#This table is constructed in 3_Calculate_parameters.R

	#Now write all data to text file for supplement
write.table(alldata, file="TableS2.txt", sep="\t")
