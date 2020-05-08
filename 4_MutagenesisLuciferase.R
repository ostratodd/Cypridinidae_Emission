require(dplyr); require(janitor)


ntalignment <- read.csv("LuciferaseTree_dNds/results/combined_codonaligned.csv",header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))

alignment <- read.csv("LuciferaseTree_dNds/results/combined_aa.csv",header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))
alignment %>% arrange(V1) -> alignment #sort by species name 
colnames(alignment)[2:ncol(alignment)] <- paste("s",seq(1,(ncol(alignment)-1)),sep=""); colnames(alignment)[1] <- "sp"
#conversion table for sites corresponding between Cypridina and present alignment
read.csv(file="LuciferaseTree_dNds/results/CypSites.csv")->cyptoalign


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
displaycodon <- function(al_sites){
    ntalignment <- read.csv("LuciferaseTree_dNds/results/combined_codonaligned.csv",header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))
	ntalignment %>% arrange(V1) -> ntalignment #sort by species name 
	
	colnames(ntalignment)[2:ncol(ntalignment)] <- paste("n",seq(1,(ncol(ntalignment)-1)),sep=""); colnames(ntalignment)[1] <- "sp"

	view <- ntalignment[,c("sp",paste("n",c((al_sites*3-2), (al_sites*3-1), (al_sites*3)), sep=""))]
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

signifmutsites <- c(38, 178, 404)
  al_signifmutsites <- cyp2aligned(signifmutsites)

	#Positive selection
al_fel <- c(43, 209) #positively selected in FEL these are ALIGNMENT numbers
  fel <- aligned2cyp(al_fel)

meme_table <- read.csv("LuciferaseTree_dNds/results/hyphy/lucclade.meme.csv",header=TRUE)
al_meme <- meme_table$Codon
  meme <- aligned2cyp(al_meme)

fel_table <- read.csv("LuciferaseTree_dNds/results/hyphy/lucclade.fel.csv",header=TRUE)

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
natural_col <- merge(alignment, tmpmerge, by='sp')
alignment_numbers <- c(38, 45, 75, 79, 87, 126, 167, 170, 178, 191, 197, 223, 258, 276, 280, 372, 375, 403, 404, 405, 406, 407, 479)
cbind(natural_col$sp,natural_col[,cyp2aligned(alignment_numbers+1)], natural_col$Lmax_Mean) -> mut_nat_col
colnames(mut_nat_col) <- colnames(mutated)

#********************Combine mutated sites and natural variants at those same sites called all_col_mut
rbind(mutated, mut_nat_col) -> all_col_mut

#*******************sites differ between Pan and Pmo
displaysites(c(38,57,110,111,171,182,187,189,224,266,273,371,478)) -> panpmo

#***********************
#Read results in csv from meme selection analysis, including positively selected sites
meme <- read.csv("LuciferaseTree_dNds/results/hyphy/lucclade.meme.csv",header=TRUE)
#Pull positively selected sites using vector of output table from meme
#Add 2 FEL sites 
c(meme_table$Codon, 43, 209) -> addfel
pos_sel <- alignment[,c("sp",paste("s",addfel, sep=""))]

#******************Now combine positively selected site data with color  data
translate <- read.csv("Raw Data/expression-kinetics/translate_lucname_decayname.csv",header=TRUE)
tmpmerge <- merge(translate, table1, by='Species')
lucNcolor <- merge(alignment, tmpmerge, by='sp')
pos_sel_col <- lucNcolor[,c("sp",paste("s",addfel, sep=""),"Lmax_Mean", "FWHM_Mean")]
remove_constant(pos_sel_col)->pos_sel_col ##Invariant sites because luc varies in spp without color data 

#******************Here combining positively selected site data with decay data
decay <- read.csv("Raw Data/expression-kinetics/decay_averages_all_for comparison_with_color.csv",header=TRUE)
translate <- read.csv("Raw Data/expression-kinetics/translate_lucname_decayname.csv",header=TRUE)
tmpmerge <- merge(decay, translate, by="PubName")
lucNkinetics <- merge(alignment, tmpmerge, by="sp")
pos_sel_decay <- lucNkinetics[,c("sp",paste("s",meme$Codon, sep=""),"lambda")]
remove_constant(pos_sel_decay)->pos_sel_decay ##Invariant sites because luc varies in spp without decay data 

#*****************Create data frame for luciferase, decay, and color for any available data
tmpmerge <- merge(translate, table1, by='Species', all=TRUE)
tmpmerge2 <- merge(tmpmerge, decay, by = 'PubName', all=TRUE)
justthesitesmaam <- displaysites(c(al_signifmutsites, al_meme, al_fel))

datamerge <- merge(justthesitesmaam, tmpmerge2, by="sp")
#remove some columns we do not need
!names(datamerge) %in% c("N", "PubName", "Species", "country", "genus", "n_lam", "Lmax_SD", "FWHM_SD", "lamSE") -> bool
datamerge[,bool] -> InterestingSites


memepvals <- lapply(meme$Episodicselectiondetected, function(x) {
               gsub("Yes_p=", "", x)
               })
c("",memepvals, "", "", "", "", "") -> mprow

data.frame(InterestingSites[-1]) -> ps_df; ps_df->ps_df_named;
rownames(ps_df_named) <- pos_sel$sp
rownames(ps_df_named) <- pos_sel$sp

round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}

#This dataframe used in Figure 1
round_df(ps_df_named, 1)->ps_df_named

#These are all the sites discussed in manuscript -- annotated
c(signifmutsites,aligned2cyp(meme$Codon), fel, "", "", "") -> cyprow
rbind(cyprow, ps_df_named) -> fullcompositetable
rownames(fullcompositetable)[1] <- "Cypridina Site Number"

rbind(fullcompositetable, c("-", "-","-","-","-","-","-","-","-","-","-","-","-", "-","-", "+", "+", "-", "-", "-", "-")) -> fullcompositetable
rownames(fullcompositetable)[17] <- "Significant Pervasive Diversifying Selection (FEL)"

rbind(fullcompositetable, c("-", "-","+","+","+","+","+","+","+","+","+","+","+", "+","+", "-", "-", "-", "-", "-", "-")) -> fullcompositetable
rownames(fullcompositetable)[18] <- "Significant Episodic Diversifying Selection (MEME)"

rbind(fullcompositetable, c("M", "M","M","-","+","+","-","-","+","-","-","-","-", "-","-", "-", "-", "-", "-", "-", "-")) -> fullcompositetable
rownames(fullcompositetable)[19] <- "Significant Correlation with Lambda Max"

rbind(fullcompositetable, c("-", "-","-","-","-","-","+","-","+","+","-","+","-", "-","-", "-", "-", "-", "-", "-", "-")) -> fullcompositetable
rownames(fullcompositetable)[20] <- "Significant Correlation with Enzymatic Decay"

fullcompositetable




