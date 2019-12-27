require(phytools)
source("plotTree.R")
luc_tree <- read.newick(file="LuciferaseTree_dNds/results/phylogenies/combined_codon.treefile")
#root(luc_tree, c("Vargula_hilgendorfii_AAA30332", "Cypridina_noctiluca_BAD08210"), resolve.root=TRUE) -> luc_tree_root
ladderize(luc_tree, right=TRUE) -> luc_tree
reroot(luc_tree, 17, resolve.root=TRUE, .5) -> luc_tree_root #root tree at midpoint of branch, which ape does not like to do
#Here is aligned amino acid file used for codon alignment. 
alignment <- read.csv("LuciferaseTree_dNds/results/combined_aa.csv",header=FALSE, stringsAsFactors=FALSE, colClasses = c("character"))
#Columns are +1 compared to meme due to sp colum Next command alters numbers by 1 to account for this and name species column sp
colnames(alignment)[2:ncol(alignment)] <- paste("s",seq(1,(ncol(alignment)-1)),sep=""); colnames(alignment)[1] <- "sp"
#Read results in csv from meme selection analysis, including positively selected sites
meme <- read.csv("LuciferaseTree_dNds/results/hyphy/lucclade.meme.csv",header=TRUE)
#Pull positively selected sites using vector of output table from meme
pos_sel <- alignment[,c("sp",paste("s",meme$Codon, sep=""))]
data.frame(pos_sel[-1]) -> ps_df; ps_df->ps_df_named;
rownames(ps_df_named) <- pos_sel$sp


