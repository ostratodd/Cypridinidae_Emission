#First remove runs where signal:noise is too low
clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
#calculate means, SDs, and counts, grouped by abbreviation
aggregate(clean[, 9:10], list(clean$abbreviation), FUN=mean) -> means
aggregate(clean[, 9:10], list(clean $abbreviation), FUN=sd) -> sds
aggregate(clean[, 9:10], list(clean $abbreviation), FUN=length) -> counts
#Now construct data frame as table for above information
cbind(subset(means, select="Group.1"), counts$sgMax, signif(means$sgMax,4), signif(sds$sgMax,2), signif(means$sgfwhm,4), signif(sds$sgfwhm,2)) -> table2
colnames(table2)[c(1,2,3,4,5,6)] <- c("Species", "N", "Lmax_Mean", "Lmax_SD", "FWHM_Mean", "FWHM_SD")
table2
	#now write to text file
write.table(table2, file = "Table2.txt", sep="\t")
#graphical table
library(gridExtra); library(grid);
quartz("Table 2")
tt3 <- ttheme_default(
  core=list(bg_params = list(fill = blues9[1:2], col=NA),
            fg_params=list(fontface=1)),
  colhead=list(fg_params=list(col="navyblue", fontface=1)),
  rowhead=list(fg_params=list(col="white")));
grid.arrange(tableGrob(table2, theme=tt3)) -> Table2

