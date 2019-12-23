#**************Functions to translate between alignment numbers for individual amino acids
aligned2cyp <- function(alignednumber) {		#This translates the present alignment number to the # in Cypridina. Input #
	#conversion table for sites corresponding between Cypridina and present alignment
	read.csv(file="LuciferaseTree_dNds/results/CypSites.csv")->cyptoalign
	cyptoalign$Cypridina_noctiluca_BAD08210[match(alignednumber,cyptoalign$Aligned)]
	
}
cyp2aligned <- function(cypnumber) {		#This translates the number in Cyprdina to the present alignment number. Input #
	#conversion table for sites corresponding between Cypridina and present alignment
	read.csv(file="LuciferaseTree_dNds/results/CypSites.csv")->cyptoalign
	cyptoalign$Aligned[match(cypnumber,cyptoalign$Cypridina_noctiluca_BAD08210)]
	
}
displaysites <- function(al_sites){
	alignment <- read.csv("LuciferaseTree_dNds/results/combined_aa.csv",header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))
	alignment %>% arrange(V1) -> alignment #sort by species name 
	
	colnames(alignment)[2:ncol(alignment)] <- paste("s",seq(1,(ncol(alignment)-1)),sep=""); colnames(alignment)[1] <- "sp"

	view <- alignment[,c("sp",paste("s",al_sites, sep=""))]
	data.frame(view[-1]) -> view_df; view_df->view_df_named;
	rownames(view_df_named) <- view$sp
	view	
}


#****************In old versions, I hand coded mutations into a csv file. These could be error prone, so did more automated approach below
#These are therefore no longer used
#For simplicity in R coding, I have 3 separate files. If I were better at R I could easily pull subsets out of 1 large file, but alas
#So, cnplus is a file that contains all the sequence data
read.csv(file="Luciferase/Cypridina only.csv") -> cn
read.csv(file="Luciferase/Cypridina mutagenesis.csv") -> cnplus
read.csv(file="Luciferase/NonCypridina only.csv") -> nc

#*****************These are important sites -- either mutated or under positive selection
	#Mutated
allmutsites <- c(38, 45, 75, 79, 87, 126, 167, 170, 178, 191, 197, 223, 258, 276, 280, 372, 375, 403, 404, 405, 406, 407, 479)
  al_allmutsites <- cyp2aligned(allmutsites)

varmutsites <- c(38, 178, 375, 404, 405) #Cypridina numbering sites with 3 states
  al_varmutsites <- cyp2aligned(varmutsites)

	#Positive selection
al_fel <- c(43, 209) #positively selected in FEL these are alignment numbers
  fel <- aligned2cyp(al_fel)

meme_table <- read.csv("LuciferaseTree_dNds/results/hyphy/lucclade.meme.csv",header=TRUE)
al_meme <- meme_table$Codon
  meme <- aligned2cyp(al_meme)


#******************Create data frame for mutations (only lmax for mutations). data frame called mutated
#csv file containing mutants and lmax data
read.csv(file="LuciferaseTree_dNds/results/Cypridina_mutations.csv", header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))->mutseqs
colnames(mutseqs)[2:ncol(mutseqs)] <- paste("c",seq(1,(ncol(mutseqs)-1)),sep="")
colnames(mutseqs)[555] <- "lmax" #last column is 555, which contains lambda max values for mutants
mutated <- data.frame(cbind(mutseqs$V1, mutseqs$c38, mutseqs$c45, mutseqs$c75, mutseqs$c79, mutseqs$c87, mutseqs$c126, mutseqs$c167, mutseqs$c170, mutseqs$c178, mutseqs$c191, mutseqs$c197, mutseqs$c223, mutseqs$c258, mutseqs$c276, mutseqs$c280, mutseqs$c372, mutseqs$c375, mutseqs$c403, mutseqs$c404, mutseqs$c405, mutseqs$c406, mutseqs$c407, mutseqs$c479, mutseqs$lmax) )
colnames(mutated) <- c("mutant", "c38", "c45", "c75", "c79","c87", "c126", "c167", "c170", "c178", "c191", "c197", "c223", "c258", "c276", "c280", "c372", "c375", "c403", "c404", "c405", "c406", "c407", "c479", "lmax")
as.numeric(as.character(mutated$lmax))->mutated$lmax

#********************Create data frame for natural variants of mutated sites.  The data fram is called mut_nat_col for mutated sites; natural; color
#combine cypridina mutations with natural luciferases for mutated sites

displaysites(al_allmutsites) -> nat_mutsites
table1 <- read.table(file="Table1.txt", sep="\t", header=TRUE) #Read again if not executed above
translate <- read.csv("Raw Data/expression-kinetics/translate_lucname_decayname.csv",header=TRUE)
tmpmerge <- merge(translate, table1, by='Species')
lucNcolor <- merge(nat_mutsites, tmpmerge, by='sp')
natural_col <- merge(dat, tmpmerge, by='sp')
cbind(natural_col$sp,natural_col[,alignment_numbers+1], natural_col$Lmax_Mean) -> mut_nat_col
colnames(mut_nat_col) <- colnames(mutated)

#********************Combine mutated sites and natural variants at those same sites called all_col_mut
rbind(mutated, mut_nat_col) -> all_col_mut
allmutationslm <- lm(lmax ~ c38 + c178 + c375 + c404 + c405, data=all_col_mut)



