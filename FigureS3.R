#############################
## plot of lamda max and decay
tmpmerge <- merge(decay, translate, by="PubName")
col_dec <- merge(table1, tmpmerge, by="Species")

plot(Lmax_Mean~lambda,data=col_dec)
cor.test(col_dec$lambda,col_dec$Lmax_Mean)
cor.test(col_dec$lambda,col_dec$FWHM_Mean)
library(ggplot2)
FigureS3 <- ggplot(data=col_dec,aes(x=lambda,y=Lmax_Mean)) + geom_point(aes(color=genus,shape=country),size=3) +
  xlab("Ave. decay constant per species") + ylab("Ave. peak emission per species")
FigureS3