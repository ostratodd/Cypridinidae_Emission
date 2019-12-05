### statistical tests for the expression assay data ###
#source publication file first 0_Publications.R
library(lme4); library(rcompanion)

##########****************************
##HEK mammalian expression statistics###

#t.tests to compare the activity of luciferases expressed in mammalian cells to the HEK cells
mam.pmo <- maxluc[maxluc$construct == "PhM" | maxluc$construct == "HEK",]
t.test(log10(mam.pmo$light)~mam.pmo$construct)
#p-value = 0.001253

mam.svu <- maxluc[maxluc$construct == "SVUluc" | maxluc$construct == "HEK",]
t.test(log10(mam.svu$light)~mam.svu$construct)
#p-value = 4.51e-09

mam.vts <- maxluc[maxluc$construct == "VtLa" | maxluc$construct == "HEK",]
t.test(log10(mam.vts$light)~mam.vts$construct)
#p-value = 4.56e-09

p.adjust(c(0.001253,4.51e-09,4.56e-09),method = "bonf")
#corrected: 3.759e-03, 1.353e-08, 1.368e-08

##########****************************
##Pichia expression statistics###

##mixed model to account for repeated measures amongst samples##
light_data_drop <- subset.data.frame(light_data,!is.na(light_data$cps_norm))
light_data_drop$Species <- factor(light_data_drop$Species,levels = c("Pichia","CNO","KHC","SVU","VTS"))
m.rand <- lme(log10(cps_norm) ~ Species*count_type, random = ~1 | sample/Replicate, data = light_data_drop)
m.rand_null <- gls(log10(cps_norm) ~ Species*count_type, data = light_data_drop)
summary(m.rand)
#these are all in comparison to Pichia, after + luciferin as the base intercept. Thus p-values indicate whether the effect is different from that standard.
#thus it seems like generally, CNO (both before and after) is different from Pichia/after
#all species have a significant interaction term, suggesting that each luciferase activity is different from Pichia/after and that changes across before/after measures

anova(m.rand,m.rand_null)
nagelkerke(m.rand,m.rand_null)
plotNormalHistogram(residuals(m.rand))
plot(fitted(m.rand), residuals(m.rand))


##Supplemental figure of the effect sizes (only works with lme4 models)
m.rand_lmer <- lmer(log10(cps_norm) ~ Species*count_type + (1|sample/Replicate),data = light_data_drop)
ggCaterpillar(ranef(m.rand_lmer)) #significant, varying random effects across samples indicate that we should account for it

##t-tests between before/after in each species to test change in light##
cno <- subset.data.frame(light_data2,light_data2$Species == "CNO")
cno <- subset.data.frame(cno,!is.na(cno$mean_cps))
t.test(log10(cno$mean_cps)~cno$count_type)
#p-value = 0.0002923

khc <- subset.data.frame(light_data2,light_data2$Species == "KHC")
t.test(log10(khc$mean_cps)~khc$count_type)
#p-value = 0.02686

vts <- subset.data.frame(light_data2,light_data2$Species == "VTS")
t.test(log10(vts$mean_cps)~vts$count_type)
#p-value = 0.005704

