#Figure 3 - Variation in Lambda-max of emission spectra
require(ggplot2)
#read in Supplemental table with all data into alldata

clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
quartz("Figure 3", 13, 3)
fig3 <- ggplot(data=clean, aes(x=abbreviation, y=sgMax, fill=abbreviation)) + geom_boxplot() + geom_jitter() 


fig3 + facet_grid(cols = vars(locality), scales = "free_x", switch = "x", space = "free")  + ylab("Lambda max (nm)") + xlab("Species") + scale_fill_manual(values=c("green", "green","#0092ff", "#0092ff", "#0092ff", "green","green", "#0092ff","#0092ff", "#0092ff","#0092ff", "green","green", "#0092ff","#0092ff", "#0092ff","#0092ff", "#0092ff", "#0092ff","#0092ff", "green","#0092ff", "#0092ff")) + theme(legend.position = "none")
