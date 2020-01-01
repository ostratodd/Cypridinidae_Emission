#Figure 3 - Variation in Lambda-max of emission spectra
require(ggplot2)
#read in Supplemental table with all data into alldata

clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
quartz("Figure 3", 13, 3)
fig3 <- qplot(abbreviation, sgMax, data=clean, geom=c("boxplot", "jitter")) + ylab("Lambda max (nm)") + xlab("Species")
	#highlight data with rectangles
	#add green blue shading w/ rectangle (geom_rect). This uses y-co-ordinates as min and max of photeros (green) and non-photeros (blue). facet_grid sorts into panels by collection locality
t2.rect1 <- data.frame (xmin=-Inf, xmax=Inf, ymin=463.97, ymax=471.14)
t1.rect1 <- data.frame (xmin=-Inf, xmax=Inf, ymin=454.27, ymax=463.42)

fig3 + geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green", alpha=0.1, inherit.aes = FALSE) + geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.07, inherit.aes = FALSE)  + facet_grid(cols = vars(locality), scales = "free_x", switch = "x", space = "free") -> Figure3
