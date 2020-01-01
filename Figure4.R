require(phytools)

InterestingSites[complete.cases(InterestingSites), ] -> full
data.frame(cbind(full$Lmax_Mean, full$lambda))->pms
row.names(pms) <- full$sp
c("Khas", "Mchi", "MspIR", "MspSVU", "Pan", "Pmor",  "PspWLU", "Vhil")-> row.names(pms)



pruned.tree<-drop.tip(luc_tree,luc_tree$tip.label[-match(full$sp, luc_tree$tip.label)])
c("Mchi", "MspIR", "MspSVU", "Pmor", "Pan", "PspWLU", "Khas", "Vhil")->pruned.tree$tip.label

phylomorphospace(pruned.tree, pms, label="radial", xlab="Lambda max (nm)", ylab="Enzymatic Decay Constant", cx=.1)


#Correlated with PIC?
#pic(pms$X1, pruned.tree)->pic.lm
#pic(pms$X2, pruned.tree)->pic.dec
#fit.pic<-lm(pic.lm~pic.dec+0)
#fit.pic
#summary(fit.pic)