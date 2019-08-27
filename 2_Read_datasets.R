
rootdir <- "~/Desktop/Color_Variation/Raw Data/"
setwd(rootdir)
#The following datafiles are written in order of collection, including some notes for each session and specimens

#For first several rounds of data collection, we used an integrating sphere. We did not change the setup
#between data collection for the first several, and so those use the same black body radiation calibration
#at first, we didn't have the black body radiation, so we did that calibration later. Since these were done
#at different times, the specific wavelengths sampled are slightly different. To arrive at the calibration, we
#interpolated data at the precise wavelengths of the sampled data. Also, for the first several data sets, we
#took a new background reading for each specimen. We found background not to change so in later data uses one
#background across each replicate collected that day.

#*****************************
#Collection date: June 8, 2016
#**************************************************************************************************
#####Roatan data. Carried back, some alive, some dried - some died in transit but were still useful
#For dried specimens in Roatan, we dried in the sun, which seems less reliable than later techniques

############GPH
workingdir <- "6-08-2016/GPH"
#background #all different backgrounds, so specified in call to graphEmission
calibration <- "../correction factor.txt"
setwd(paste0(rootdir,workingdir))
graphEmission("GPH_1_05sec_x10x1sp_live_sig.asc", "GPH_1_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_1;
graphEmission("GPH_2_05sec_x10x1sp_live_sig.asc", "GPH_2_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_2;
graphEmission("GPH_3_05sec_x10x1sp_live_sig.asc", "GPH_3_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_3;
graphEmission("GPH_4_05sec_x10x1sp_live_sig.asc", "GPH_4_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_4;
graphEmission("GPH_5_05sec_x10x1sp_live_sig.asc", "GPH_5_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_5;
graphEmission("GPH_6_05sec_x10x1sp_live_sig.asc", "GPH_6_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_6;
graphEmission("GPH_7_05sec_x10x1sp_live_sig.asc", "GPH_7_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_7;
graphEmission("GPH_8_05sec_x10x1sp_live_sig.asc", "GPH_8_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> GPH_8;

##############WLU
workingdir <- "6-08-2016/WLU"
#background #all different backgrounds, so specified in call to graphEmission
calibration <- "../correction factor.txt"
setwd(paste0(rootdir,workingdir))
graphEmission("WLU_1_05sec_x5x1sp_sig.asc", "WLU_1_05sec_x5x1sp_bgd.asc", calibration,5,TRUE) -> WLU_1;
graphEmission("WLU_2_05sec_x5x1sp_sig.asc", "WLU_2_05sec_x5x1sp_bgd.asc", calibration,5,TRUE) -> WLU_2;
graphEmission("WLU_3_05sec_x5x1sp_sig.asc", "WLU_3_05sec_x5x1sp_bgd.asc", calibration,5,TRUE) -> WLU_3;
graphEmission("WLU_4_05sec_x5x1sp_sig.asc", "WLU_4_05sec_x5x1sp_bgd.asc", calibration,5,TRUE) -> WLU_4;
graphEmission("WLU_5_05sec_x5x1sp_sig.asc", "WLU_5_05sec_x5x1sp_bgd.asc", calibration,5,TRUE) -> WLU_5;
graphEmission("WLU_6_05sec_x5x1sp_sig.asc", "WLU_6_05sec_x5x1sp_bgd.asc", calibration,5,TRUE) -> WLU_6;
#graphEmission("WLU_7_05sec_x10x1sp_dry_sig.asc", "WLU_7_05sec_x10x1sp_dry_bgd.asc", calibration,10,TRUE) -> WLU_7; NO signal here
graphEmission("WLU_8_05sec_x10x3sp_dry_sig.asc", "WLU_8_05sec_x10x3sp_dry_bgd.asc", calibration,10,TRUE) -> WLU_8;
graphEmission("WLU_10_05sec_x10x3sp_dry_sig.asc", "WLU_10_05sec_x10x3sp_dry_bgd.asc", calibration,10,TRUE) -> WLU_10;

##############LSD
workingdir <- "6-08-2016/LSD"
setwd(paste0(rootdir,workingdir))
#background #all different backgrounds, so specified in call to graphEmission
calibration <- "../correction factor.txt"
graphEmission("LSD1.txt", "background.txt", calibration,10,TRUE) -> LSD_1;
graphEmission("LSD2.txt", "background.txt", calibration,10,TRUE) -> LSD_2;

##############IR
workingdir <- "6-08-2016/IR"
setwd(paste0(rootdir,workingdir))
#background #all different backgrounds, so specified in call to graphEmission
calibration <- "../correction factor.txt"
graphEmission("IR_1_05sec_x5x1sp_dead_sig.asc", "IR_1_05sec_x5x1sp_dead_bgd.asc", calibration,10,TRUE) -> IR_1;
graphEmission("IR_2_05sec_x5x1sp_dead_sig.asc", "IR_2_05sec_x5x1sp_dead_bgd.asc", calibration,10,TRUE) -> IR_2;
graphEmission("IR_3_05sec_x5x1sp_live_sig.asc", "IR_3_05sec_x5x1sp_live_bgd.asc", calibration,10,TRUE) -> IR_3;
graphEmission("IR_4_05sec_x5x1sp_dead_sig.asc", "IR_4_05sec_x5x1sp_dead_bgd.asc", calibration,10,TRUE) -> IR_4;
graphEmission("IR_5_05sec_x5x1sp_dead_sig.asc", "IR_5_05sec_x5x1sp_dead_bgd.asc", calibration,10,TRUE) -> IR_5;
graphEmission("IR_6_05sec_x5x1sp_dead_sig.asc", "IR_6_05sec_x5x1sp_dead_bgd.asc", calibration,10,TRUE) -> IR_6;

##############RD
workingdir <- "6-08-2016/RD"
setwd(paste0(rootdir,workingdir))
#background #all different backgrounds, so specified in call to graphEmission
calibration <- "../correction factor.txt"
graphEmission("RD_1_05sec_x5x1sp_dead_sig.asc", "RD_1_05sec_x5x1sp_dead_bgd.asc", calibration,5,TRUE) -> RD_1;
graphEmission("RD_2_05sec_x5x1sp_dead_sig.asc", "RD_2_05sec_x5x1sp_dead_bgd.asc", calibration,5,TRUE) -> RD_2;
graphEmission("RD_3_05sec_x5x1sp_dead_sig.asc", "RD_3_05sec_x5x1sp_dead_bgd.asc", calibration,5,TRUE) -> RD_3;
graphEmission("RD_4_05sec_x5x1sp_dead_sig.asc", "RD_4_05sec_x5x1sp_dead_bgd.asc", calibration,5,TRUE) -> RD_4;
graphEmission("RD_5_05sec_x5x1sp_dead_sig.asc", "RD_5_05sec_x5x1sp_dead_bgd.asc", calibration,5,TRUE) -> RD_5;

###############DU
workingdir <- "6-08-2016/DU"
setwd(paste0(rootdir,workingdir))
#background #all different backgrounds, so specified in call to graphEmission
calibration <- "../correction factor.txt"

graphEmission("DU_1_05sec_x10x1sp_live_sig.asc", "DU_1_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> DU_1;
graphEmission("DU_2_05sec_x10x1sp_live_sig.asc", "DU_2_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> DU_2;
graphEmission("DU_3_05sec_x10x1sp_live_sig.asc", "DU_3_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> DU_3;
graphEmission("DU_4_05sec_x10x1sp_live_sig.asc", "DU_4_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> DU_4;
graphEmission("DU_5_05sec_x10x1sp_live_sig.asc", "DU_5_05sec_x10x1sp_live_bgd.asc", calibration,10,TRUE) -> DU_5;


#********************************
#Collection date: July 14, 2016
#Vargula tsujii from Catalina Island
##############Vts -- Vargula tsujii from Catalina Island, CA; USA
workingdir <- "7-14-2016/Vtsu"
setwd(paste0(rootdir,workingdir))
#background #all different backgrounds, so specified in call to graphEmission
calibration <- "../correction factor.txt"
graphEmission("Vtsu_fem1_sig.asc", "Vtsu_fem1_bgd.asc", calibration,10,TRUE) -> Vtsu_1;
graphEmission("Vtsu_fem2_sig.asc", "Vtsu_fem2_bgd.asc", calibration,10,TRUE) -> Vtsu_2;
graphEmission("Vtsu_fem3_sig.asc", "Vtsu_fem3_bgd.asc", calibration,10,TRUE) -> Vtsu_3;
graphEmission("Vtsu_male1_sig.asc", "Vtsu_male1_bgd.asc", calibration,5,TRUE) -> Vtsu_4;
graphEmission("Vtsu_male2_sig.asc", "Vtsu_male2_bgd.asc", calibration,5,TRUE) -> Vtsu_5;
graphEmission("Vtsu_male3_sig.asc", "Vtsu_male3_bgd.asc", calibration,5,TRUE) -> Vtsu_6;
graphEmission("Vtsu_male4_sig.asc", "Vtsu_male4_bgd.asc", calibration,10,TRUE) -> Vtsu_7;
graphEmission("Vtsu_male5_sig.asc", "Vtsu_male5_bgd.asc", calibration,10,TRUE) -> Vtsu_8;


#********************************
#Collection date: August 10, 2016
#Puerto Rico specimens 8-10-2016*******************************************************************************

##############WCU
workingdir <- "8-10-2016/WCU"
setwd(paste0(rootdir,workingdir))
background <- "background.txt"
calibration <- "../correction factor.txt"

graphEmission("WCU1.txt", background, calibration,10,TRUE) -> WCU_1;
graphEmission("WCU2.txt", "background.txt", calibration,10,TRUE) -> WCU_2;
graphEmission("WCU3.txt", "background.txt", calibration,10,TRUE) -> WCU_3;
#graphEmission("WCU4.txt", "background.txt", calibration,10,TRUE) -> WCU4; #WCU3 and WCU4 are the same data files
graphEmission("WCU5.txt", "background.txt", calibration,10,TRUE) -> WCU_5;

##############SPU -- A Kornickeria from Puerto Rico
workingdir <- "8-10-2016/SPU"
setwd(paste0(rootdir,workingdir))
background <- "background.asc"
calibration <- "../correction factor.txt"
graphEmission("SPU1_v2.txt", background, calibration,10,TRUE) -> SPU_1;
graphEmission("SPU2.asc.txt", background, calibration,10,TRUE) -> SPU_2;
graphEmission("SPU3.asc", background, calibration,10,TRUE) -> SPU_3;
graphEmission("SPU4.asc", background, calibration,10,TRUE) -> SPU_4;
background <- "background2.txt"
graphEmission("SPU5.txt", background, calibration,10,TRUE) -> SPU_5;
graphEmission("SPU6.txt", background, calibration,10,TRUE) -> SPU_6;
graphEmission("SPU7.txt", background, calibration,10,TRUE) -> SPU_7;

#**********************************
#Data Collection: September 7, 2016
#********************************************************************
#Ordered dried Vargula hilgendorfii from Carolina Biological. This is
#to compare with spectrum taken by Yasuo Mitani in Japan, and to published
#spectra of Vargula hilgendorfii. Also, used these dried samples later to 
#serve as a sort of standard.
#Vhilg
workingdir <- "9-07-2016/Vhilg"
setwd(paste0(rootdir,workingdir))
background <- "bgd.asc"
calibration <- "../correction factor.txt"

graphEmission("Vhilg1.asc", background, calibration,10,TRUE) -> Vhil090716_1;
graphEmission("Vhilg2.asc", background, calibration,10,TRUE) -> Vhil090716_2;
graphEmission("Vhilg3.asc", background, calibration,10,TRUE) -> Vhil090716_3;
graphEmission("Vhilg4.asc", background, calibration,10,TRUE) -> Vhil090716_4;
graphEmission("Vhilg5_DeleteLastAcq.asc", background, calibration,9,TRUE) -> Vhil090716_5;
graphEmission("Vhilg6.asc", background, calibration,10,TRUE) -> Vhil090716_6;


#**********************************
#Data Collection: September 20, 2016
#********************************************************************
#Ordered dried Vargula hilgendorfii from Carolina Biological. 
#Vhilg
workingdir <- "9-20-2016/Vhilg"
setwd(paste0(rootdir,workingdir))
background <- "background_bgd.asc"
calibration <- "../correction factor.txt"

graphEmission("Vhilg1_3indiv.asc", background, calibration,10,TRUE) -> Vhil092016_1;
graphEmission("Vhilg_2_3indiv.asc", background, calibration,10,TRUE) -> Vhil092016_2;
graphEmission("Vhilg3_3indiv.asc", background, calibration,10,TRUE) -> Vhil092016_3;




#**********************************
#Data Collection: October 5, 2016
#********************************************************************
#From this point on, we removed integrating sphere to read from cuvette for higher sensitivity
##############Photeros annecohenae sent by Karl from near Belize City. Different population than normal Belize Pa, so added "karl" to end
workingdir <- "10-05-2016/Pann"
setwd(paste0(rootdir,workingdir))
#background is different for each specimen
calibration <- "../correction factor nosphere.txt"

graphEmission("PhotAnnKarl1_sig.asc", "PhotAnnKarl1_bgd.asc", calibration,10,TRUE) -> Pann_1;
graphEmission("PhotAnnKarl2_sig.asc", "PhotAnnKarl2_bgd.asc", calibration,10,TRUE) -> Pann_2;
graphEmission("PhotAnnKarl3_sig.asc", "PhotAnnKarl3_bgd.asc", calibration,10,TRUE) -> Pann_3;

#Vhilg as standard. 
workingdir <- "10-05-2016/Vhilg"
setwd(paste0(rootdir,workingdir))
#background is different for each specimen
calibration <- "../correction factor nosphere.txt"
graphEmission("VhilgendorfiiOcto5_sig.asc", "VhilgendorfiiOcto5_bgd.asc", calibration,10,TRUE) -> Vhil10052016_1;


#****************************
#Data Collection: January 23, 2017
#********************************************************************
#####Belize data 1-23-17 sent by Gretchen Gerrish to UCSB. Dried.

background <- "../1-23-17background.asc";
calibration <- "../correction factor.txt";
	#Khastingsi-------------------------------------------------------------
workingdir <- "1-23-2017/Khastingsi"
setwd(paste0(rootdir,workingdir))
graphEmission(  "KhastingsiDry.asc", background, calibration,20,TRUE) -> Khas_1;
graphEmission("KhastingsiLive1.asc", background, calibration,20,TRUE) -> Khas_2;
graphEmission("KhastingsiLive2.asc", background, calibration,20,TRUE) -> Khas_3;

	#MSH-------------------------------------------------------------
workingdir <- "1-23-2017/MSH"
setwd(paste0(rootdir,workingdir))
graphEmission(  "MSHlive1.asc", background, calibration,20,TRUE) -> MSH_1;
graphEmission(  "MSHlive2.asc", background, calibration,20,TRUE) -> MSH_2;
graphEmission(  "MSHlive3.asc", background, calibration,20,TRUE) -> MSH_3;
graphEmission(  "MSHlive4.asc", background, calibration,20,TRUE) -> MSH_4;
graphEmission(  "MSHlive5.asc", background, calibration,20,TRUE) -> MSH_5;

	#Pmor-------------------------------------------------------------
workingdir <- "1-23-2017/Pmorini"
setwd(paste0(rootdir,workingdir))
graphEmission(  "PmoriniLive1.asc", background, calibration,20,TRUE) -> Pmor_1;
graphEmission(  "PmoriniLive2.asc", background, calibration,20,TRUE) -> Pmor_2;
graphEmission(  "PmoriniLive3.asc", background, calibration,20,TRUE) -> Pmor_3;
graphEmission(  "PmoriniLive4.asc", background, calibration,20,TRUE) -> Pmor_4;
graphEmission(  "Pmorini5remax.txt", background, calibration,17,TRUE) -> Pmor_5; #first 3 maxed sensor, removed

	#SVD-------------------------------------------------------------
workingdir <- "1-23-2017/SVD"
setwd(paste0(rootdir,workingdir))
graphEmission(  "SVDLive1.asc", background, calibration,20,TRUE) -> SVD_1;
graphEmission(  "SVDLive2.asc", background, calibration,20,TRUE) -> SVD_2;
graphEmission(  "SVDLive3.asc", background, calibration,20,TRUE) -> SVD_3;
graphEmission(  "SVDLive4.asc", background, calibration,20,TRUE) -> SVD_4;
graphEmission(  "SVDLive5.asc", background, calibration,20,TRUE) -> SVD_5;

	#SVU-------------------------------------------------------------
workingdir <- "1-23-2017/SVU"
setwd(paste0(rootdir,workingdir))
graphEmission(  "SVULive1.asc", background, calibration,20,TRUE) -> SVU_1;
graphEmission(  "SVULive2.asc", background, calibration,20,TRUE) -> SVU_2;
graphEmission(  "SVULive3.asc", background, calibration,20,TRUE) -> SVU_3;
graphEmission(  "SVULive4.asc", background, calibration,20,TRUE) -> SVU_4;
graphEmission(  "SVUremax5.txt", background, calibration,16,TRUE) -> SVU_5;	#first 4 maxed sensor, removed

#Vargula hilgendorfii control-------------------------------------------------------------
workingdir <- "1-23-2017/Vhil_standard"
setwd(paste0(rootdir,workingdir))
graphEmission(  "VhilgendorfiiDry.asc", background, calibration,20,TRUE) -> Vhil012317_1;
graphEmission(  "VhilgendorfiiDry2.asc", background, calibration,20,TRUE) -> Vhil012317_2;
graphEmission(  "VhilgendorfiiDry3.asc", background, calibration,20,TRUE) -> Vhil012317_3;

#*****************************
#Collection Date: May 30, 2017
#********************************************************************
#5-30-2017
#Dried species from Panama May 2017
# collected at Bocas del Toro, and dried in oven. Carried with silica to keep dry
workingdir <- "5-30-2017/EGD"
background <- "../baselinePanamaSeriesMay26.asc";
calibration <- "../panama species correction factor.txt"

#EGD-------------------------------------------------------------
setwd(paste0(rootdir,workingdir))
graphEmission("EGD_dried_male1_removeMaxed.txt", background, calibration,13,TRUE) -> EGD_1; #First 2 maxed, so deleted
#graphEmission("EGD_dried_male2.asc", background, calibration, 15, TRUE) -> EGD_2; #no light
graphEmission("EGD_dried_male3.asc", background, calibration,15,TRUE) -> EGD_3;
graphEmission("EGD_dried_male4.asc", background, calibration,15,TRUE) -> EGD_4; #check 2 strange outliers
graphEmission("EGD_dried_male5.asc", background, calibration,15,TRUE) -> EGD_5;

#Contragula------------------------------------------------------
workingdir <- "5-30-2017/CONT"
setwd(paste0(rootdir,workingdir))
graphEmission("CONT_Panama_dried_male1_removeMaxed.txt", background, calibration,12,TRUE) -> Cont_bocas_1; #removed first 3 scans maxed sensor
graphEmission("CONT_Panama_dried_male2.asc", background, calibration,15,TRUE) -> Cont_bocas_2;
graphEmission("CONT_Panama_dried_male3.asc", background, calibration,15,TRUE) -> Cont_bocas_3;

#LLL-------------------------------------------------------------
workingdir <- "5-30-2017/LLL"
setwd(paste0(rootdir,workingdir))
graphEmission("LLL_dried_male1.asc", background, calibration,15,TRUE) -> LLL_1;
graphEmission("LLL_dried_male2.asc", background, calibration,15,TRUE) -> LLL_2;
graphEmission("LLL_dried_male3.asc", background, calibration,15,TRUE) -> LLL_3;
graphEmission("LLL_dried_male4.asc", background, calibration,15,TRUE) -> LLL_4;
graphEmission("LLL_dried_male5.asc", background, calibration,15,TRUE) -> LLL_5;

#MFU-------------------------------------------------------------
workingdir <- "5-30-2017/MFU"
setwd(paste0(rootdir,workingdir))
graphEmission("MFU_male1_remax.txt", background, calibration,14,TRUE) -> MFU_1; #First sample is maxed, so deleted
graphEmission("MFU_dried_male2.asc", background, calibration,15,TRUE) -> MFU_2;
graphEmission("MFU_dried_male3.asc", background, calibration,15,TRUE) -> MFU_3;
graphEmission("MFU_male4_remax.txt", background, calibration,12,TRUE) -> MFU_4; #First 3samples maxed, so deleted
graphEmission("MFU_dried_male5.asc", background, calibration,15,TRUE) -> MFU_5;

#SFM-------------------------------------------------------------
workingdir <- "5-30-2017/SFM"
setwd(paste0(rootdir,workingdir))
graphEmission("SFM_dried_male1.asc", background, calibration,15,TRUE) -> SFM_1;
graphEmission("SFM_dried_male2.asc", background, calibration,15,TRUE) -> SFM_2;
graphEmission("SFM_dried_male3.asc", background, calibration,15,TRUE) -> SFM_3;
graphEmission("SFM_dried_male4.asc", background, calibration,15,TRUE) -> SFM_4;
graphEmission("SFM_dried_male5.asc", background, calibration,15,TRUE) -> SFM_5;

#SMU-------------------------------------------------------------
workingdir <- "5-30-2017/SMU"
setwd(paste0(rootdir,workingdir))
graphEmission("SMU_dried_male1.asc", background, calibration,15,TRUE) -> SMU_1;
graphEmission("SMU_dried_male2.asc", background, calibration,15,TRUE) -> SMU_2;
graphEmission("SMU_dried_male3.asc", background, calibration,15,TRUE) -> SMU_3;
graphEmission("SMU_male4_remax.txt", background, calibration,11,TRUE) -> SMU_4; #last 4 maxed to deleted
graphEmission("SMU_dried_male5.asc", background, calibration,15,TRUE) -> SMU_5;

#*******************************************************************************
#Non-ucsb data and published data
#Emily Ellis and Yasuo Mitani set up integrating sphere and spectroradiometer to record emission data in Japan
#This data is already calculated in photons for each wavelength, so no need to run graphEmission to subtract background
workingdir <- "non-ucsb/Vhil"
setwd(paste0(rootdir,workingdir))
raw <- read.csv(file="Vh.csv");  
Vhun <- data.frame(wavelength = raw$Wavelength.nm., sum=raw$Var.1.21sec.);
Vhil_Japan <- data.frame(wavelength = Vhun$wavelength, sum=Vhun$sum/max(Vhun$sum, na.rm=T));

#Published spectrum for Vargula hilgendorfii
raw <- read.table(file="Vh from Tsuji paper Tess scan.txt");
Vhun <- data.frame(wavelength = raw$V1, sum=raw$V2);
Vhil_tsuji <- data.frame(wavelength = Vhun$wavelength, sum=Vhun$sum/max(Vhun$sum, na.rm=T));

#Published spectrum for Cypridina noctiluca in patent by Ohmiya
workingdir <- "non-ucsb/Cnoc"
setwd(paste0(rootdir,workingdir))
raw <- read.table(file="noctiluca WT datathief.txt");
Cnoc_ohmiya <- data.frame(wavelength = raw$V1, sum=raw$V2);
#################Pg -- Photeros gramminicola from Huvard 1993, digitized - already smoothed.
workingdir <- "non-ucsb/Pgra"
setwd(paste0(rootdir,workingdir))
raw <- read.table(file="Pg.txt");
Pgra_huvard <- data.frame(wavelength = raw$V1, sum=raw$V2);
