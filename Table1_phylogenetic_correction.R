### phylogenetic ANOVA methods for the color paper ###
library(geiger); library(nlme); library(phytools)

#read in color data and positively selected sites from 4_MutagenesisLuciferase.R first
head(pos_sel_col)

#read in gene tree
luc_tree <- read.newick(file="LuciferaseTree_dNds/results/phylogenies/combined_codon.treefile")

##analysis for all the color site data including mutants from Table 3.R
head(all_col_mut)
all_ct <- read.tree("LuciferaseTree_dNds/results/phylogenies/color_data_all_mutants_tree.txt")
#match tips and datasets
col_seq <- all_ct$tip.label
nomat_col <- col_seq[!(col_seq %in% all_col_mut$mutant)]
pgls_col <- drop.tip(all_ct,nomat_col) #create var-covar matrix with BM model
dt2 <- all_col_mut
row.names(dt2) <- dt2$mutant
dt2 = dt2[match(pgls_col$tip.label, row.names(dt2)),] #data and tree must be in same order
bm.model_2 = corBrownian(1,pgls_col)

#drop factors not signifcant from initial ANOVA model c372 & c405, and interactions except for c375*c404
col_mdl_RE <- gls(lmax ~ c178 + c280 + c375 * c404,data=dt2,correlation = bm.model_2,method="REML",control = list(singular.ok = TRUE))
col_mdl_ML <- gls(lmax ~ c178 + c280 + c375 * c404,data=dt2,correlation = bm.model_2,method="ML",control = list(singular.ok = TRUE))
col_mdl_no <- gls(lmax ~ c178 + c280 + c375 * c404,data=dt2,method="ML",control = list(singular.ok = TRUE))
AICc(col_mdl_ML); AICc(col_mdl_no)
summary(col_mdl_no)
f <- anova(col_mdl_no)

#this analysis is for Table4.R that restricts to the last three color sites
tree_seq <- luc_tree$tip.label
nomat_tr <- tree_seq[!(tree_seq %in% pos_sel_col$sp)]
pgls_tree <- drop.tip(luc_tree,nomat_tr) #create var-covar matrix with BM model
dt3 <- pos_sel_col
row.names(dt3) <- dt3$sp
dt3 = dt3[match(pgls_tree$tip.label, row.names(dt3)),] #data and tree must be in same order
bm.model_3 = corBrownian(1,pgls_tree)

lmax_bm_est <- gls(Lmax_Mean ~ s102 + s142 + s189,correlation = bm.model_3, data = dt3,method="REML")
lmax_bm <- gls(Lmax_Mean ~ s102 + s142 + s189,correlation = bm.model_3, data = dt3,method="ML")
lmax_no <- gls(Lmax_Mean ~ s102 + s142 + s189, data = dt3,method="ML")
AIC(lmax_bm); AIC(lmax_no)
anova(lmax_no, type = "marginal")

#this analysis is for Table5.R that restricts to the three kinetics sites
dec_seq <- luc_tree$tip.label
nomat <- dec_seq[!(dec_seq %in% pos_sel_decay$sp)]
pgls_dec <- drop.tip(luc_tree,nomat) #create var-covar matrix with BM model
dt4 <- pos_sel_decay
row.names(dt4) <- dt4$sp
dt4 = dt4[match(pgls_dec$tip.label, row.names(dt4)),] #data and tree must be in same order
bm.model_4 = corBrownian(1,pgls_dec)

decay_bm_RE <- gls(lambda ~ s41 + s189 + s261, correlation = bm.model_4, data = dt4,method="REML",control = list(singular.ok = TRUE))
decay_bm_ML <- gls(lambda ~ s41 + s189 + s261, correlation = bm.model_4, data = dt4,method="ML",control = list(singular.ok = TRUE))
decay_bm_no <- gls(lambda ~ s41 + s189 + s261, data = dt4,method="ML",control = list(singular.ok = TRUE))
AICc(decay_bm_ML); AICc(decay_bm_no) #prefers no correction
anova(decay_bm_no)
#removing site s41 as it is non-significant in a GLS framework
decay_bm_RE2 <- gls(lambda ~ s189 + s261, correlation = bm.model_4, data = dt4,method="REML",control = list(singular.ok = TRUE))
decay_bm_ML2 <- gls(lambda ~ s189 + s261, correlation = bm.model_4, data = dt4,method="ML",control = list(singular.ok = TRUE))
decay_bm_no2 <- gls(lambda ~ s189 + s261, data = dt4,method="ML",control = list(singular.ok = TRUE))
AICc(decay_bm_ML2); AICc(decay_bm_no2)
#models are equivalent with or without phylogenetic correction by AICc
anova(decay_bm_ML2,test = "marginal")
anova(decay_bm_no2,test="marginal")
