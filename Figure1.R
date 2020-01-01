require(phytools)
source("plotTree.R")
luc_tree <- read.newick(file="LuciferaseTree_dNds/results/phylogenies/combined_codon.treefile")
#root(luc_tree, c("Vargula_hilgendorfii_AAA30332", "Cypridina_noctiluca_BAD08210"), resolve.root=TRUE) -> luc_tree_root
ladderize(luc_tree, right=TRUE) -> luc_tree
reroot(luc_tree, 17, resolve.root=TRUE, .5) -> luc_tree_root #root tree at midpoint of branch, which ape does not like to do




#executing commands below in separate block yield tree of desired proportion
quartz("Figure 1", 10, 4)
plotTree_mod(tree=luc_tree_root,ancestral.reconstruction=F,tip.labels=TRUE, lwd=1, infoFile=ps_df_named, treeWidth=7,infoWidth=4, infoCols=c("Lmax_Mean", "s41", "s93", "s102", "s142", "s160", "s177", "s189", "s261", "s285", "s291", "s320", "s389", "s477", "lambda"))
