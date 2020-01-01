#Figure S2 - Variation in FWHM of emission spectra
require(ggplot2)
#Read in supplemental table into variable alldata

clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")

quartz("Figure S2", 13, 3)
FigureS2 <- qplot(abbreviation, sgfwhm, data=clean, geom=c("boxplot", "jitter")) + ylab("FWHM (nm)") + xlab("Species")
FigureS2 + facet_grid(cols = vars(locality), scales = "free_x", switch = "x", space = "free") -> FigureS2
FigureS2

