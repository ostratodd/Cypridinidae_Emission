#Figure 3 - Variation in Lambda-max of emission spectra
require(ggplot2)
#read in Supplemental table with all data into alldata

clean <- subset(alldata, error < 0.02 | replicate=="Vhil_tsuji" | replicate=="Cnoc_ohmiya" | replicate=="Pgra_huvard")
quartz("Figure 3", 13, 3)
fig3 <- ggplot(data=clean, aes(x= reorder(label, -sgMax, fun = median,), y=sgMax, fill=abbreviation)) + geom_boxplot() + geom_jitter()

clean$label[clean$abbreviation == "P_EGD"] <- "P. sp. EGD (Pa)"
clean$label[clean$abbreviation == "P_SFM"] <- "P. sp. SFM (Pa)"
clean$label[clean$abbreviation == "M_CON"] <- "C. sp. CONT (Pa)"
clean$label[clean$abbreviation == "M_LLL"] <- "M. sp. LLL (Pa)"
clean$label[clean$abbreviation == "M_MFU"] <- "M. sp. MFU (Pa)"
clean$label[clean$abbreviation == "P_ann"] <- "P. annecohenae (Bz)"
clean$label[clean$abbreviation == "P_mor"] <- "P. morini (Bz)"
clean$label[clean$abbreviation == "M_SMU"] <- "M. sp. SMU (Pa)"
clean$label[clean$abbreviation == "K_has"] <- "K. hastingsi (Bz)"
clean$label[clean$abbreviation == "M_chi"] <- "M. chicoi (Bz)"
clean$label[clean$abbreviation == "M_SVU"] <- "M. sp. SVU (Bz)"
clean$label[clean$abbreviation == "P_GPH"] <- "P. sp. GPH (Ro)"
clean$label[clean$abbreviation == "P_WLU"] <- "P. sp. WLU (Ro)"
clean$label[clean$abbreviation == "M_LSD"] <- "M. sp. LSD (PR)"
clean$label[clean$abbreviation == "M_IR"] <- "P. sp. IR (Ro)"
clean$label[clean$abbreviation == "M_RD"] <- "P. sp. RD (Ro)"
clean$label[clean$abbreviation == "M_DU"] <- "P. sp. DU (Ro)"
clean$label[clean$abbreviation == "V_tsu"] <- "V. tsujii (US)"
clean$label[clean$abbreviation == "K_SPU"] <- "K. sp. SPU (PR)"
clean$label[clean$abbreviation == "V_hil"] <- "V. hilgendorfii (Jp)"
clean$label[clean$abbreviation == "V_hil_p"] <- "V. hilgendorfii * (Jp)"
clean$label[clean$abbreviation == "P_gra_p"] <- "P. graminicola * (Pa)"
clean$label[clean$abbreviation == "C_noc_p"] <- "C. noctiluca * (Jp)"

fig3 + facet_grid(cols = vars(genus), scales = "free_x", switch = "x", space = "free")  + ylab("Lambda max (nm)") + xlab("Species") +
  scale_fill_manual(values=c("#00a9ff", "#00a9ff","#007bff", "#007bff", "#007bff", "#00a9ff","#00a9ff", "#007bff","#007bff", "#007bff","#007bff", "#00a9ff","#00a9ff", "#007bff","#007bff", "#007bff","#007bff", "#007bff", "#007bff","#007bff", "#00a9ff","#007bff", "#007bff")) +
  theme(legend.position = "none",axis.text.x = element_text(angle = -60,hjust = 0,)) -> Figure3

