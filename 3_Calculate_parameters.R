#First execute functions

#Next read in data


#Here, calculate parameters and add metadata
#***********************************EGD
locality <- "Panama"
genus <- "Photeros"
species<-"Photeros_EGD"
abbreviation <- "P_EGD"
sex<-"male"
specloc <- "ucsb"
preservation <- "dried"

fullrow(abbreviation, locality, genus,  species, "EGD1", EGD_1, sex, preservation, specloc) -> EGD1f
fullrow(abbreviation, locality, genus,  species, "EGD3", EGD_3, sex, preservation, specloc) -> EGD3f
fullrow(abbreviation, locality, genus,  species, "EGD4", EGD_4, sex, preservation, specloc) -> EGD4f
fullrow(abbreviation, locality, genus,  species, "EGD5", EGD_5, sex, preservation, specloc) -> EGD5f

#Now bind all samples together into one large data frame
rbind(EGD1f, EGD3f, EGD4f, EGD5f) -> allEGD

#***********************************EGD
genus <- "Contragula"
species<-"contragula"
abbreviation<-"CONT"
sex<-"male"
specloc <- "ucsb"

fullrow(abbreviation, locality, genus,  species, "cont1", Cont_bocas_1, sex, preservation, specloc) -> cont1f
fullrow(abbreviation, locality, genus,  species, "cont2", Cont_bocas_2, sex, preservation, specloc) -> cont2f
fullrow(abbreviation, locality, genus,  species, "cont3", Cont_bocas_3, sex, preservation, specloc) -> cont3f

#Now bind all samples together into one large data frame
rbind(cont1f, cont2f, cont3f) -> allcont

#***********************************LLL
genus <- "Maristella"
species<-"LLL"
abbreviation<-"LLL"
sex<-"male"
specloc <- "ucsb"

fullrow(abbreviation, locality, genus,  species, "LLL1", LLL_1, sex, preservation, specloc) -> LLL1f
fullrow(abbreviation, locality, genus,  species, "LLL2", LLL_2, sex, preservation, specloc) -> LLL2f
fullrow(abbreviation, locality, genus,  species, "LLL3", LLL_3, sex, preservation, specloc) -> LLL3f
fullrow(abbreviation, locality, genus,  species, "LLL4", LLL_4, sex, preservation, specloc) -> LLL4f
fullrow(abbreviation, locality, genus,  species, "LLL5", LLL_5, sex, preservation, specloc) -> LLL5f

rbind(LLL1f, LLL2f, LLL3f, LLL4f, LLL5f) -> allLLL


#***********************************MFU
genus <- "Maristella"
species<-"MFU"
abbreviation<-"MFU"
sex<-"male"
specloc <- "ucsb"

fullrow(abbreviation, locality, genus,  species, "MFU1", MFU_1, sex, preservation, specloc) -> MFU1f
fullrow(abbreviation, locality, genus,  species, "MFU2", MFU_2, sex, preservation, specloc) -> MFU2f
fullrow(abbreviation, locality, genus,  species, "MFU3", MFU_3, sex, preservation, specloc) -> MFU3f
fullrow(abbreviation, locality, genus,  species, "MFU4", MFU_4, sex, preservation, specloc) -> MFU4f
fullrow(abbreviation, locality, genus,  species, "MFU5", MFU_5, sex, preservation, specloc) -> MFU5f

rbind(MFU1f, MFU2f, MFU3f, MFU4f, MFU5f) -> allMFU

#***********************************SFM
genus <- "Photeros"
species<-"Photeros_SFM"
abbreviation<-"SFM"
sex<-"male"
specloc <- "ucsb"

fullrow(abbreviation, locality, genus,  species, "SFM1", SFM_1, sex, preservation, specloc) -> SFM1f
fullrow(abbreviation, locality, genus,  species, "SFM2", SFM_2, sex, preservation, specloc) -> SFM2f
fullrow(abbreviation, locality, genus,  species, "SFM3", SFM_3, sex, preservation, specloc) -> SFM3f
fullrow(abbreviation, locality, genus,  species, "SFM4", SFM_4, sex, preservation, specloc) -> SFM4f
fullrow(abbreviation, locality, genus,  species, "SFM5", SFM_5, sex, preservation, specloc) -> SFM5f

rbind(SFM1f, SFM2f, SFM3f, SFM4f, SFM5f) -> allSFM

#***********************************SFM
locality <- "Belize"
genus <- "Maristella"
species<-"Maristella_SMU"
abbreviation<-"SMU"
sex<-"male"
specloc <- "ucsb"

fullrow(abbreviation, locality, genus,  species, "SMU1", SMU_1, sex, preservation, specloc) -> SMU1f
fullrow(abbreviation, locality, genus,  species, "SMU2", SMU_2, sex, preservation, specloc) -> SMU2f
fullrow(abbreviation, locality, genus,  species, "SMU3", SMU_3, sex, preservation, specloc) -> SMU3f
fullrow(abbreviation, locality, genus,  species, "SMU4", SMU_4, sex, preservation, specloc) -> SMU4f
fullrow(abbreviation, locality, genus,  species, "SMU5", SMU_5, sex, preservation, specloc) -> SMU5f

rbind(SMU1f, SMU2f, SMU3f, SMU4f, SMU5f) -> allSMU


#***********************************Kornickeria hastingsi
genus <- "Kornickeria"
species<-"Kornickeria_hastingsi"
abbreviation<-"K_has"
sex<-"male"
specloc <- "ucsb"
preservation <- "dried"
fullrow(abbreviation, locality, genus,  species, "Khas1", Khas_1, sex, preservation, specloc) -> Khas1f
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "Khas2", Khas_2, sex, preservation, specloc) -> Khas2f
fullrow(abbreviation, locality, genus,  species, "Khas3", Khas_3, sex, preservation, specloc) -> Khas3f

rbind(Khas1f,Khas2f,Khas3f) -> allKhas

#***********************************Maristella chicoi
genus <- "Maristella"
species<-"Maristella_chicoi"
abbreviation<-"M_chi"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "MSH1", MSH_1, sex, preservation, specloc) -> MSH1f
fullrow(abbreviation, locality, genus,  species, "MSH2", MSH_2, sex, preservation, specloc) -> MSH2f
fullrow(abbreviation, locality, genus,  species, "MSH3", MSH_3, sex, preservation, specloc) -> MSH3f
fullrow(abbreviation, locality, genus,  species, "MSH4", MSH_4, sex, preservation, specloc) -> MSH4f
fullrow(abbreviation, locality, genus,  species, "MSH5", MSH_5, sex, preservation, specloc) -> MSH5f

rbind(MSH1f,MSH2f,MSH3f,MSH4f,MSH5f) -> allMSH

#***********************************Photeros morini
genus <- "Photeros"
species<-"Photeros_morini"
abbreviation<-"P_mor"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "Pmor1", Pmor_1, sex, preservation, specloc) -> Pmor1f
fullrow(abbreviation, locality, genus,  species, "Pmor2", Pmor_2, sex, preservation, specloc) -> Pmor2f
fullrow(abbreviation, locality, genus,  species, "Pmor3", Pmor_3, sex, preservation, specloc) -> Pmor3f
fullrow(abbreviation, locality, genus,  species, "Pmor4", Pmor_4, sex, preservation, specloc) -> Pmor4f
fullrow(abbreviation, locality, genus,  species, "Pmor5", Pmor_5, sex, preservation, specloc) -> Pmor5f

rbind(Pmor1f,Pmor2f,Pmor3f,Pmor4f,Pmor5f) -> allPmor

#***********************************SVD
genus <- "Maristella"
species<-"Maristella_SVD"
abbreviation<-"M_SFU"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "SVD1", SVD_1, sex, preservation, specloc) -> SVD1f
fullrow(abbreviation, locality, genus,  species, "SVD2", SVD_2, sex, preservation, specloc) -> SVD2f
fullrow(abbreviation, locality, genus,  species, "SVD3", SVD_3, sex, preservation, specloc) -> SVD3f
fullrow(abbreviation, locality, genus,  species, "SVD4", SVD_4, sex, preservation, specloc) -> SVD4f
fullrow(abbreviation, locality, genus,  species, "SVD5", SVD_5, sex, preservation, specloc) -> SVD5f

rbind(SVD1f,SVD2f,SVD3f,SVD4f,SVD5f) -> allSVD

#***********************************SVD
genus <- "Maristella"
species<-"Maristella_SVU"
abbreviation<-"M_SVU"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "SVU1", SVU_1, sex, preservation, specloc) -> SVU1f
fullrow(abbreviation, locality, genus,  species, "SVU2", SVU_2, sex, preservation, specloc) -> SVU2f
fullrow(abbreviation, locality, genus,  species, "SVU3", SVU_3, sex, preservation, specloc) -> SVU3f
fullrow(abbreviation, locality, genus,  species, "SVU4", SVU_4, sex, preservation, specloc) -> SVU4f
fullrow(abbreviation, locality, genus,  species, "SVU5", SVU_5, sex, preservation, specloc) -> SVU5f

rbind(SVU1f,SVU2f,SVU3f,SVU4f,SVU5f) -> allSVU

#***********************************Vargula hilgendorfii
genus <- "Vargula"
locality <- "Purchased"
species<-"Vargula_hilgendorfii"
abbreviation<-"V_hil"
sex<-"unknown"
specloc <- "ucsb"
preservation <- "dried"
fullrow(abbreviation, locality, genus,  species, "Vhil0123171", Vhil012317_1, sex, preservation, specloc) -> Vhil0123171f
fullrow(abbreviation, locality, genus,  species, "Vhil0123172", Vhil012317_2, sex, preservation, specloc) -> Vhil0123172f
fullrow(abbreviation, locality, genus,  species, "Vhil0123173", Vhil012317_3, sex, preservation, specloc) -> Vhil0123173f

rbind(Vhil0123171f, Vhil0123172f, Vhil0123173f) -> allVhil012317

species<-"Vargula_hilgendorfii"
abbreviation<-"V_hil"
sex<-"unknown"
specloc <- "ucsb"
preservation <- "dried"
fullrow(abbreviation, locality, genus,  species, "Vhil0907161", Vhil090716_1, sex, preservation, specloc) -> Vhil0907161f
fullrow(abbreviation, locality, genus,  species, "Vhil0907162", Vhil090716_2, sex, preservation, specloc) -> Vhil0907162f
fullrow(abbreviation, locality, genus,  species, "Vhil0907163", Vhil090716_3, sex, preservation, specloc) -> Vhil0907163f
fullrow(abbreviation, locality, genus,  species, "Vhil0907164", Vhil090716_4, sex, preservation, specloc) -> Vhil0907164f
fullrow(abbreviation, locality, genus,  species, "Vhil0907165", Vhil090716_5, sex, preservation, specloc) -> Vhil0907165f
fullrow(abbreviation, locality, genus,  species, "Vhil0907166", Vhil090716_6, sex, preservation, specloc) -> Vhil0907166f

	rbind(Vhil0907161f, Vhil0907162f, Vhil0907163f, Vhil0907164f, Vhil0907165f, Vhil0907166f) -> allVhil090716

fullrow(abbreviation, locality, genus,  species, "Vhil0920161", Vhil092016_1, sex, preservation, specloc) -> Vhil0920161f
fullrow(abbreviation, locality, genus,  species, "Vhil0920162", Vhil092016_2, sex, preservation, specloc) -> Vhil0920162f
fullrow(abbreviation, locality, genus,  species, "Vhil0920163", Vhil092016_3, sex, preservation, specloc) -> Vhil0920163f
	rbind(Vhil0920161f, Vhil0920162f, Vhil0920163f) -> allVhil092016

	#one trial from 10-05 
fullrow(abbreviation, locality, genus,  species, "Vhil100520161", Vhil10052016_1, sex, preservation, specloc) -> Vhil100520161f
	rbind(Vhil100520161f) -> allVhil10052016
specloc <- "Japan"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "Vhil_Japan", Vhil_Japan, sex, preservation, specloc) -> VhilJapan1f
	rbind(VhilJapan1f) -> allVhilJapan

#***********************************GPH
locality <- "Roatan"
genus <- "Photeros"
species<-"Photeros_GPH"
abbreviation<-"P_GPH"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"

fullrow(abbreviation, locality, genus,  species, "GPH1", GPH_1, sex, preservation, specloc) -> GPH1f
fullrow(abbreviation, locality, genus,  species, "GPH2", GPH_2, sex, preservation, specloc) -> GPH2f
fullrow(abbreviation, locality, genus,  species, "GPH3", GPH_3, sex, preservation, specloc) -> GPH3f
fullrow(abbreviation, locality, genus,  species, "GPH4", GPH_4, sex, preservation, specloc) -> GPH4f
fullrow(abbreviation, locality, genus,  species, "GPH5", GPH_5, sex, preservation, specloc) -> GPH5f
fullrow(abbreviation, locality, genus,  species, "GPH6", GPH_6, sex, preservation, specloc) -> GPH6f
fullrow(abbreviation, locality, genus,  species, "GPH7", GPH_7, sex, preservation, specloc) -> GPH7f
fullrow(abbreviation, locality, genus,  species, "GPH8", GPH_8, sex, preservation, specloc) -> GPH8f

rbind(GPH1f, GPH2f, GPH3f, GPH4f, GPH5f, GPH6f, GPH7f, GPH8f) -> allGPH
genus <- "Photeros"
species<-"Photeros_WLU"
abbreviation<-"P_WLU"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "WLU1", WLU_1, sex, preservation, specloc) -> WLU1f
fullrow(abbreviation, locality, genus,  species, "WLU2", WLU_2, sex, preservation, specloc) -> WLU2f
fullrow(abbreviation, locality, genus,  species, "WLU3", WLU_3, sex, preservation, specloc) -> WLU3f
fullrow(abbreviation, locality, genus,  species, "WLU4", WLU_4, sex, preservation, specloc) -> WLU4f
fullrow(abbreviation, locality, genus,  species, "WLU5", WLU_5, sex, preservation, specloc) -> WLU5f
fullrow(abbreviation, locality, genus,  species, "WLU6", WLU_6, sex, preservation, specloc) -> WLU6f
preservation <- "dried"
fullrow(abbreviation, locality, genus,  species, "WLU7", WLU_8, sex, preservation, specloc) -> WLU8f
fullrow(abbreviation, locality, genus,  species, "WLU8", WLU_10, sex, preservation, specloc) -> WLU10f

rbind(WLU1f, WLU2f, WLU3f, WLU4f, WLU5f, WLU6f, WLU8f, WLU10f) -> allWLU
genus <- "Maristella"
species<-"LSD"
abbreviation<-"LSD"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "LSD1", LSD_1, sex, preservation, specloc) -> LSD1f
fullrow(abbreviation, locality, genus,  species, "LSD2", LSD_2, sex, preservation, specloc) -> LSD2f

rbind(LSD1f, LSD2f) -> allLSD

#***********************************IR
genus <- "Maristella"
species<-"IR"
abbreviation<-"IR"
sex<-"male"
specloc <- "ucsb"
preservation <- "dead"

fullrow(abbreviation, locality, genus,  species, "IR1", IR_1, sex, preservation, specloc) -> IR1f
fullrow(abbreviation, locality, genus,  species, "IR2", IR_2, sex, preservation, specloc) -> IR2f
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "IR3", IR_3, sex, preservation, specloc) -> IR3f
preservation <- "dead"
fullrow(abbreviation, locality, genus,  species, "IR4", IR_4, sex, preservation, specloc) -> IR4f
fullrow(abbreviation, locality, genus,  species, "IR5", IR_5, sex, preservation, specloc) -> IR5f
fullrow(abbreviation, locality, genus,  species, "IR6", IR_6, sex, preservation, specloc) -> IR6f

rbind(IR1f, IR2f, IR3f, IR4f, IR5f, IR6f) -> allIR

genus <- "Maristella"
species<-"RD"
abbreviation<-"RD"
sex<-"male"
specloc <- "ucsb"
preservation <- "dead"
fullrow(abbreviation, locality, genus,  species, "RD1", RD_1, sex, preservation, specloc) -> RD1f
fullrow(abbreviation, locality, genus,  species, "RD2", RD_2, sex, preservation, specloc) -> RD2f
fullrow(abbreviation, locality, genus,  species, "RD3", RD_3, sex, preservation, specloc) -> RD3f
fullrow(abbreviation, locality, genus,  species, "RD4", RD_4, sex, preservation, specloc) -> RD4f
fullrow(abbreviation, locality, genus,  species, "RD5", RD_5, sex, preservation, specloc) -> RD5f

rbind(RD1f, RD2f, RD3f, RD4f, RD5f) -> allRD

genus <- "Maristella"
species<-"DU"
abbreviation<-"DU"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "DU1", DU_1, sex, preservation, specloc) -> DU1f
fullrow(abbreviation, locality, genus,  species, "DU2", DU_2, sex, preservation, specloc) -> DU2f
fullrow(abbreviation, locality, genus,  species, "DU3", DU_3, sex, preservation, specloc) -> DU3f
fullrow(abbreviation, locality, genus,  species, "DU4", DU_4, sex, preservation, specloc) -> DU4f
fullrow(abbreviation, locality, genus,  species, "DU5", DU_5, sex, preservation, specloc) -> DU5f

rbind(DU1f, DU2f, DU3f, DU4f, DU5f) -> allDU

#***********************************Vargula tsujii
locality <- "California_USA"
genus <- "Fred"
species<-"Vargula_tsujii"
abbreviation<-"V_tsu"
sex<-"female"
specloc <- "ucsb"
preservation <- "live"

fullrow(abbreviation, locality, genus,  species, "Vtsu1", Vtsu_1, sex, preservation, specloc) -> Vtsu1f
fullrow(abbreviation, locality, genus,  species, "Vtsu2", Vtsu_2, sex, preservation, specloc) -> Vtsu2f
fullrow(abbreviation, locality, genus,  species, "Vtsu3", Vtsu_3, sex, preservation, specloc) -> Vtsu3f
sex<-"male"
fullrow(abbreviation, locality, genus,  species, "Vtsu4", Vtsu_4, sex, preservation, specloc) -> Vtsu4f
fullrow(abbreviation, locality, genus,  species, "Vtsu5", Vtsu_5, sex, preservation, specloc) -> Vtsu5f
fullrow(abbreviation, locality, genus,  species, "Vtsu6", Vtsu_6, sex, preservation, specloc) -> Vtsu6f
fullrow(abbreviation, locality, genus,  species, "Vtsu7", Vtsu_7, sex, preservation, specloc) -> Vtsu7f
fullrow(abbreviation, locality, genus,  species, "Vtsu8", Vtsu_8, sex, preservation, specloc) -> Vtsu8f

rbind(Vtsu1f, Vtsu2f, Vtsu3f, Vtsu4f, Vtsu5f, Vtsu6f, Vtsu7f, Vtsu8f) -> allVtsu


locality <- "Puerto Rico_USA"
genus <- "Kornickeria"
species<-"Kornickeria_WCU"
abbreviation<-"K_WCU"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "WCU1", WCU_1, sex, preservation, specloc) -> WCU1f
fullrow(abbreviation, locality, genus,  species, "WCU2", WCU_2, sex, preservation, specloc) -> WCU2f
fullrow(abbreviation, locality, genus,  species, "WCU3", WCU_3, sex, preservation, specloc) -> WCU3f
fullrow(abbreviation, locality, genus,  species, "WCU5", WCU_5, sex, preservation, specloc) -> WCU5f

rbind(WCU1f, WCU2f, WCU3f, WCU5f) -> allWCU

#***********************************SPU
genus <- "Kornickeria"
species<-"Kornickeria_SPU"
abbreviation<-"K_SPU"
sex<-"male"
specloc <- "ucsb"
preservation <- "live"

fullrow(abbreviation, locality, genus,  species, "SPU1", SPU_1, sex, preservation, specloc) -> SPU1f
fullrow(abbreviation, locality, genus,  species, "SPU2", SPU_2, sex, preservation, specloc) -> SPU2f
fullrow(abbreviation, locality, genus,  species, "SPU3", SPU_3, sex, preservation, specloc) -> SPU3f
fullrow(abbreviation, locality, genus,  species, "SPU4", SPU_4, sex, preservation, specloc) -> SPU4f
fullrow(abbreviation, locality, genus,  species, "SPU5", SPU_5, sex, preservation, specloc) -> SPU5f
fullrow(abbreviation, locality, genus,  species, "SPU6", SPU_6, sex, preservation, specloc) -> SPU6f
fullrow(abbreviation, locality, genus,  species, "SPU7", SPU_7, sex, preservation, specloc) -> SPU7f

rbind(SPU1f, SPU2f, SPU3f, SPU4f, SPU5f, SPU6f, SPU7f) -> allSPU

#***********************************Photeros annecohenae
locality <- "Belize"
genus <- "Photeros"
species<-"Photeros_annecohenae"
abbreviation<-"P_ann"
sex<-"unknown"
specloc <- "ucsb"
preservation <- "live"
fullrow(abbreviation, locality, genus,  species, "Pann1", Pann_1, sex, preservation, specloc) -> Pann1f
fullrow(abbreviation, locality, genus,  species, "Pann2", Pann_2, sex, preservation, specloc) -> Pann2f
fullrow(abbreviation, locality, genus,  species, "Pann3", Pann_3, sex, preservation, specloc) -> Pann3f

rbind(Pann1f,Pann2f,Pann3f) -> allPann

#***********************************Previously published data
locality <- "Published"
genus <- "Vargula"
species<-"Vargula_hilgendorfii"
abbreviation<-"V_hil_p"
sex<-"unknown"
specloc <- "published"
preservation <- "unknown"

fullrow(abbreviation, locality, genus,  species, "Vhil_tsuji", Vhil_tsuji, sex, preservation, specloc) -> Vhiltsujif
locality <- "Published"
genus <- "Cypridina"
species<-"Cypridina_noctiluca"
abbreviation<-"C_noc_p"
sex<-"unknown"
specloc <- "published"
preservation <- "unknown"

fullrow(abbreviation, locality, genus,  species, "Cnoc_ohmiya", Cnoc_ohmiya, sex, preservation, specloc) -> Cnocohmiyaf
locality <- "Published"
genus <- "Photeros"
species<-"Photeros_gramminicola"
abbreviation<-"P_gra_p"
sex<-"unknown"
specloc <- "published"
preservation <- "unknown"

fullrow(abbreviation, locality, genus,  species, "Pgra_huvard", Pgra_huvard, sex, preservation, specloc) -> Pgrahuvardf



#****************************************Put all data into a single data frame, using one row bind command to bind them all
#onering
rbind(allEGD,allcont,allLLL, allMFU, allSFM, allSMU, allKhas, allMSH, allPmor, allSVD, allSVU, allVhil012317, allVhil090716, allVhil10052016, allVhil092016, allVhilJapan, allGPH, allWLU, allLSD, allIR, allRD, allDU, allVtsu, allWCU, allSPU, allPann, Vhiltsujif, Cnocohmiyaf, Pgrahuvardf) -> alldata
alldata

