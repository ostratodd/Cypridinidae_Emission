
#Plot dn:ds p-values across gene data from hyphy is read from csv files in 4_ file
#For FEL plot, should remove sites with low dn:ds, since p-value could represent high or low dn:ds

plot(data=all_fel_table, (1-p.value)~Site, ylab="(1-p value)", xlab='Aligned Codon Site', ylim=c(0.85,1), cex=0.3, col="dark grey", type="h", pch=19)
par(new=TRUE)
plot(data=all_meme_table, (1-p.value)~Site, ylab="",xlab='Aligned Codon Site', ylim=c(0.85,1), cex=0.3, col="blue", type="h", pch=19)
