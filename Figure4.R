require(phytools)
#****************Function for plotting the ancestral states
plotASRsite <- function(site, title) {
	x <- alignment[, site+1]
	data.frame(x, row.names=alignment$sp) -> x

	cursite <- setNames(x[,1], rownames(x))
	fitER<-ace(cursite,luc_tree_root,model="ER",type="discrete")
	fitER
	fitER$lik.anc
	cols<-setNames(c("red","blue", "green", "orange", "yellow", "black"),levels(cursite))
	quartz(title, 8, 5)
	plotTree(luc_tree_root,fsize=.6,ftype="i",lwd=1)

	nodelabels(node=1:luc_tree_root$Nnode+Ntip(luc_tree_root),
    		pie=fitER$lik.anc,piecol=cols,cex=0.4)
	tiplabels(pie=to.matrix(cursite[luc_tree_root$tip.label],
    		levels(cursite)),piecol=cols,cex=0.3)	
} 
#******************* are all the ancestral state reconstructions, which were moved by hand to phylophenospace
ladderize(luc_tree, right=TRUE) -> luc_tree
reroot(luc_tree, 17, resolve.root=TRUE, .5) -> luc_tree_root #root tree at midpoint of branch, which ape does not like to do

plotASRsite(207, "Site 178 - Neutrally evolving color site")
plotASRsite(cyp2aligned(280), "Site 280 - Neutrally evolving color site")
plotASRsite(cyp2aligned(375), "Site 375 - Neutrally evolving color site")


plotASRsite(435, "Site 404 - Color site under negative selection")
plotASRsite(102, "Site 74 - Correlated with color and under positive selection")



plotASRsite(189, "Site 160 - Correlated with color and decay and under positive selection")

plotASRsite(41, "Site 19 - Correlated with decay and under positive selection")
plotASRsite(142, "Site 114 - Correlated with decay and under positive selection")
plotASRsite(cyp2aligned(232), "Site 232 - Correlated with decay and under positive selection")


#******************Phylophenospace starts here
InterestingSites[complete.cases(InterestingSites), ] -> full
data.frame(cbind(full$Lmax_Mean, full$lambda))->pms
row.names(pms) <- full$sp
c("Khas", "Mchi", "MspIR", "MspSVU", "Pan", "Pmor",  "PspWLU", "Vhil")-> row.names(pms)



pruned.tree<-drop.tip(luc_tree,luc_tree$tip.label[-match(full$sp, luc_tree$tip.label)])
c("Mchi", "MspIR", "MspSVU", "Pmor", "Pan", "PspWLU", "Khas", "Vhil")->pruned.tree$tip.label

quartz("Figure 4 - Phylophenospace (Final modified in Illustrator)")
phylomorphospace(pruned.tree, pms, label="radial", xlab="Lambda max (nm)", ylab="Enzymatic Decay Constant", cx=.1)


#*****************Shows no correlation between phenotypes, even when using PICs
#Correlated with PIC?
#pic(pms$X1, pruned.tree)->pic.lm
#pic(pms$X2, pruned.tree)->pic.dec
#fit.pic<-lm(pic.lm~pic.dec+0)
#fit.pic
#summary(fit.pic)



