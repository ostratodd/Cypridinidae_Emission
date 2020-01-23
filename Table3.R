################
#Table 3 - ANOVA results for 5 ‘mutagenesis sites’ using all data from mutagenesis study and available species’ luciferases

#aligned versus cypridina sites
#64  71 103 107 115 154 196 199 207 220 226 252 287 307 311 403 406 434 435 436 437 438 512
#38  45  75  79  87 126 167 170 178 191 197 223 258 276 280 372 375 403 404 405 406 407 479

anova(lm(lmax ~ c38 * c178 * c375 * c404 * c405, data=all_col_mut)) -> table3
write.table(table3, file = "Table3.txt", sep="\t")
table3_mutagenesis

#Model selection does not choose c38
anova(lm(lmax ~ c178 * c280 * c372 * c375 * c404 * c405, data=all_col_mut)) -> table3





##Below uses model selection, instead of using sites with most variants
#This takes a while to run, so commented. Uncomment to re-run
#library(gridExtra); library(grid); require(janitor);
#library(MuMIn);
#
#options(na.action = "na.fail")
#amutationslm <- lm(lmax ~ c38+ c45+ c75+ c79+ c87+ c126+ c167+ c170+ c178+ c191+ c197+ c223+ c258+ c276+ c280+ c372+ c375+ c403+ c404+ c405+ c406+ c407+ c479, data=all_col_mut)
#mixnmatch_amut <- dredge(amutationslm,rank = "AIC",m.lim = c(0,6)) 
#head(mixnmatch_amut, 20)


#Results of model selection
#
#Global model call: lm(formula = lmax ~ c38 + c45 + c75 + c79 + c87 + c126 + c167 + 
#    c170 + c178 + c191 + c197 + c223 + c258 + c276 + c280 + c372 + 
#    c375 + c403 + c404 + c405 + c406 + c407 + c479, data = all_col_mut)
#---
#Model selection table 
#        (Intrc) c126 c167 c170 c178 c191 c197 c223 c276 c280 c372 c375 c403 c404 c405 c406 c407 c45 c479 c75 c79 df  logLik   AIC delta weight
#51977     454.6                   +                   +    +         +         +    +                            15 -78.291 186.6  0.00  0.052
#52489     454.6                   +                   +         +    +         +    +                            15 -78.291 186.6  0.00  0.052
#51721     454.2                   +                        +         +         +    +                            14 -79.461 186.9  0.34  0.044
#52233     454.2                   +                             +    +         +    +                            14 -79.461 186.9  0.34  0.044
#52745     454.2                   +                        +    +    +         +    +                            14 -79.461 186.9  0.34  0.044
#51785     455.9                   +              +         +         +         +    +                            15 -78.632 187.3  0.68  0.037
#52297     455.9                   +              +              +    +         +    +                            15 -78.632 187.3  0.68  0.037
#1100297   455.9                   +                        +         +         +    +                      +     15 -78.632 187.3  0.68  0.037
#1100809   455.9                   +                             +    +         +    +                      +     15 -78.632 187.3  0.68  0.037
#51722     455.9    +              +                        +         +         +    +                            15 -78.632 187.3  0.68  0.037
#52234     455.9    +              +                             +    +         +    +                            15 -78.632 187.3  0.68  0.037
#2148873   454.4                   +                        +         +         +    +                          + 15 -78.870 187.7  1.16  0.029
#2149385   454.4                   +                             +    +         +    +                          + 15 -78.870 187.7  1.16  0.029
#59913     454.6                   +                        +         +    +    +    +                            15 -78.879 187.8  1.17  0.029
#60425     454.6                   +                             +    +    +    +    +                            15 -78.879 187.8  1.17  0.029
#
#
#51737     454.2                   +    +                   +         +         +    +                            15 -79.340 188.7  2.10  0.018
#
