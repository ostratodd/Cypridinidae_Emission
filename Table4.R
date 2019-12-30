colorlm <- lm(Lmax_Mean ~ s41 + s93 + s102 + s142 + s160 + s177 + s189 + s261 + s285 + s291 + s320 + s389 + s477 , data=pos_sel_col)
      
mixnmatch_col <- dredge(	colorlm,rank = "AIC",m.lim = c(0,3)) #3 seems like the maximum terms we can fit safely
head(mixnmatch_col, 7)

#Global model call: lm(formula = Lmax_Mean ~ s41 + s93 + s102 + s142 + s160 + s177 + 
#    s189 + s261 + s285 + s291 + s320 + s389 + s477, data = pos_sel_col)
#---
#Model selection table 
#     (Intrc) s102 s142 s189 s261 s285 s320 s93 df logLik  AIC delta weight
#51     458.8         +    +    +                9 10.990 -4.0  0.00  0.167
#113    458.8              +    +    +           9 10.990 -4.0  0.00  0.167
#305    458.8              +    +         +      9 10.990 -4.0  0.00  0.167
#82     458.8    +         +         +           9 10.990 -4.0  0.00  0.167
#274    458.8    +         +              +      9 10.990 -4.0  0.00  0.167
#20     458.8    +    +    +                     9 10.990 -4.0  0.00  0.167


#4131   458.0         +         +             +  7 -7.838 29.7 33.66  0.000

lmax_anova <- lm(Lmax_Mean ~ s102 + s142 + s189 + s261,data=pos_sel_col) # rep'd 3 or more x in best models
anova(lmax_anova)	#102, 142, 189
lmax_anova2 <- lm(Lmax_Mean ~ s102 + s142 + s189 + s261 + s285 + s320, data=pos_sel_col) # rep'd 2 or more x in best models
anova(lmax_anova2)	#102, 142, 189 significant
displaysites(c(102,142,189)) -> positivelyselectedThatPredictLmax
positivelyselectedThatPredictLmax
aligned2cyp(c(102,142,189))

table4 <- anova(lmax_anova)



fwhmlm <- lm(FWHM_Mean ~ s41 + s93 + s102 + s142 + s160 + s177 + s189 + s261 + s285 + s291 + s320 + s389 + s477 , data=pos_sel_col)
