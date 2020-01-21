require(phytools)
source("plotTree.R")
luc_tree <- read.newick(file="LuciferaseTree_dNds/results/phylogenies/combined_codon.treefile")

luc_tree_boot <- read.newick(file="LuciferaseTree_dNds/results/phylogenies/combined_codon.contree")

ladderize(luc_tree, right=TRUE) -> luc_tree
reroot(luc_tree, 17, resolve.root=TRUE, .5) -> luc_tree_root #root tree at midpoint of branch, which ape does not like to do
ladderize(luc_tree_boot, right=TRUE) -> luc_tree_boot
reroot(luc_tree_boot, 17, resolve.root=TRUE, .5) -> luc_tree_boot_root #root tree at midpoint of branch, which ape does not like to do

#plotTree function to put table next to  tree does not seem to print node label
quartz("Bootstraps 1", 10, 4)
plot(luc_tree_boot_root, show.node.label=TRUE)

#executing commands below in separate block yield tree of desired proportion
quartz("Figure 1", 10, 4)
plotTree_mod(tree=luc_tree_boot_root,ancestral.reconstruction=F,tip.labels=TRUE,  lwd=1, infoFile=ps_df_named, treeWidth=7,infoWidth=4, infoCols=c("Lmax_Mean", "s41", "s93", "s102", "s142", "s160", "s177", "s189", "s261", "s285", "s291", "s320", "s389", "s477", "lambda"))

#Not all interesting sites currently in ps_df_named
#s64 s207 s435 s41 s93 s102 s142 s160 s177 s189 s261 s285 s291 s320 s389 s477 s43 s209 Lmax_Mean FWHM_Mean lambda