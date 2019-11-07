rootdir <- paste0(maindir,"Luciferase/")	#maindir is set in calling function for figures 0_Publication_figures.R. rootdir is where rawdata originates for the current script
setwd(rootdir)

#For simplicity in R coding, I have 3 separate files. If I were better at R I could easily pull subsets out of 1 large file, but alas
#So, cnplus is a file that contains all the sequence data
read.csv(file="Cypridina only.csv") -> cn
read.csv(file="Cypridina mutagenesis.csv") -> cnplus
read.csv(file="NonCypridina only.csv") -> nc





