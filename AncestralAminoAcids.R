require(phytools)

luc_tree <- read.newick(file="LuciferaseTree_dNds/results/phylogenies/combined_codon.treefile")
ladderize(luc_tree, right=TRUE) -> luc_tree
reroot(luc_tree, 17, resolve.root=TRUE, .5) -> luc_tree_root #root tree at midpoint of branch, which ape does not like to do

par(mfrow=c(2,2))
plotASRsite(207, "Site 207 - Neutrally evolving color site")
plotASRsite(404, "Site 404 - Color site under negative selection")

plotASRsite(102, "Site 102 - Correlated with color and decay and under positive selection")
plotASRsite(189, "Site 189 - Correlated with color and decay and under positive selection")


plotASRsite <- function(site, title) {
	x <- alignment[, site+1]
	data.frame(x, row.names=alignment$sp) -> x

	cursite <- setNames(x[,1], rownames(x))
	fitER<-ace(cursite,luc_tree_root,model="ER",type="discrete")
	fitER
	fitER$lik.anc

	quartz(title, 8, 5)
	plotTree(luc_tree_root,fsize=.6,ftype="i",lwd=1)

	nodelabels(node=1:luc_tree_root$Nnode+Ntip(luc_tree_root),
    		pie=fitER$lik.anc,piecol=cols,cex=0.4)
	tiplabels(pie=to.matrix(cursite[luc_tree_root$tip.label],
    		levels(cursite)),piecol=cols,cex=0.3)

	
	
	
} 
