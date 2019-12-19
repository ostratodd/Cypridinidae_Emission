rootdir <- paste0(maindir,"Luciferase/")	#maindir is set in calling function for figures 0_Publication_figures.R. rootdir is where rawdata originates for the current script
setwd(rootdir)

#For simplicity in R coding, I have 3 separate files. If I were better at R I could easily pull subsets out of 1 large file, but alas
#So, cnplus is a file that contains all the sequence data
read.csv(file="Cypridina only.csv") -> cn
read.csv(file="Cypridina mutagenesis.csv") -> cnplus
read.csv(file="NonCypridina only.csv") -> nc

#csv file containing mutants and lmax data
read.csv(file="../LuciferaseTree_dNds/results/Cypridina_mutations.csv", header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))->mutseqs
colnames(mutseqs)[2:ncol(mutseqs)] <- paste("c",seq(1,(ncol(mutseqs)-1)),sep="")
colnames(mutseqs)[555] <- "lmax" #last column is 555, which contains lambda max values for mutants
allmutsites <- c(38, 45, 75, 79, 87, 126, 167, 170, 178, 191, 197, 223, 258, 276, 280, 372, 375, 403, 404, 405, 406, 407, 479)
mutsites <- c(38, 178, 375, 404, 405) #Cypridina numbering sites with 3 states

mutated <- data.frame(cbind(mutseqs$V1, mutseqs$c38, mutseqs$c45, mutseqs$c75, mutseqs$c79, mutseqs$c87, mutseqs$c126, mutseqs$c167, mutseqs$c170, mutseqs$c178, mutseqs$c191, mutseqs$c197, mutseqs$c223, mutseqs$c258, mutseqs$c276, mutseqs$c280, mutseqs$c372, mutseqs$c375, mutseqs$c403, mutseqs$c404, mutseqs$c405, mutseqs$c406, mutseqs$c407, mutseqs$c479, mutseqs$lmax) )
colnames(mutated) <- c("mutant", "c38", "c45", "c75", "c79","c87", "c126", "c167", "c170", "c178", "c191", "c197", "c223", "c258", "c276", "c280", "c372", "c375", "c403", "c404", "c405", "c406", "c407", "c479", "lmax")


#conversion table for sites corresponding between Cypridina and present alignment
read.csv(file="../LuciferaseTree_dNds/results/CypSites.csv")->cyptoalign

maindir <- "~/Documents/GitHub/Cypridinidae_EmissionSpectra/"
setwd(maindir)
