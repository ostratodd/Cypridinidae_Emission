### statistical tests for the expression assay data ###
#source publication file first 0_Publications.R
library(lme4); library(rcompanion)

##########****************************
##HEK mammalian expression statistics###




##########****************************
##Pichia expression statistics###

##mixed model to account for repeated measures amongst samples##
light_data_drop <- subset.data.frame(light_data,!is.na(light_data$cps_norm))
light_data_drop$Species <- factor(light_data_drop$Species,levels = c("Pichia","CNO","KHC","SVU","VTS"))
m.rand <- lme(log10(cps_norm) ~ Species*count_type, random = ~1 | sample/Replicate, data = light_data_drop)
m.rand_null <- gls(log10(cps_norm) ~ Species*count_type, data = light_data_drop)
summary(m.rand)
anova(m.rand,m.rand_null)
nagelkerke(m.rand,m.rand_null)
plotNormalHistogram(residuals(m.rand))
plot(fitted(m.rand), residuals(m.rand))

##Supplemental figure of the effect sizes (only works with lme4 models)
m.rand_lmer <- lmer(log10(cps_norm) ~ Species*count_type + (1|sample/Replicate),data = light_data_drop)
ggCaterpillar(ranef(m.rand_lmer)) #variable random effects across samples indicate that we should account for it

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

p.adjust(c(0.0002923,0.02686,0.005704),method = "bonf")

