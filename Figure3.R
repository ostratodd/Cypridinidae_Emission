#Figure 3 - Variation in Lambda-max of emission spectra
require(ggplot2)
#read in Supplemental table with all data into alldata

clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
quartz("Figure 3", 13, 3)
fig3 <- ggplot(data=clean, aes(x=abbreviation, y=sgMax, fill=abbreviation)) + geom_boxplot() + geom_jitter() 


fig3 + facet_grid(cols = vars(genus), scales = "free_x", switch = "x", space = "free")  + ylab("Lambda max (nm)") + xlab("Species") + scale_fill_manual(values=c("#00a9ff", "#00a9ff","#007bff", "#007bff", "#007bff", "#00a9ff","#00a9ff", "#007bff","#007bff", "#007bff","#007bff", "#00a9ff","#00a9ff", "#007bff","#007bff", "#007bff","#007bff", "#007bff", "#007bff","#007bff", "#00a9ff","#007bff", "#007bff")) + theme(legend.position = "none")
