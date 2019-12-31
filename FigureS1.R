require(ape)
require(phytools)
quartz("Figure S1")
read.tree(file="./LuciferaseTree_dNds/results/phylogenies/all_aa.treefile") -> lucplus_tree
midpoint.root(lucplus_tree) -> lucplus_tree_r
ladderize(lucplus_tree_r)-> lprl
plot(lprl, show.tip.label=TRUE, cex=.8, x.lim=2)
