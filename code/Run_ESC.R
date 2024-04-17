
# date: 17-04-24
# author: VB
# desc: script to run ESC for x native tree spp, for baseline & three climate scenarios

### working dirs ----

#wd <- "~/Documents/Woodland-Trust/Data-Analysis/SoWT2_Spp-uncertainty" # Mac path
wd <- "C:/Users/vbu/OneDrive - the Woodland Trust/Data-analysis/SoWT2_Spp-uncertainty" # WT path
dirData <-  paste0(wd,"/data-raw/")
dirScratch <-  paste0(wd,"/data-scratch/")
dirOut <- paste0(wd,"/data-out/")


### libs ----

library(tidyverse)
library(raster)
library(rgdal)


### ESC parameters and function ----

# amend raster options - adjust memory etc 
rasterOptions(maxmemory=2e+09, chunksize=1e+06,overwrite=T,tmptime=1.1)

# read in parameters
parameters <- read.csv(paste0(dirData,"ESC/esccoeffs_crafty.csv"))
Species <- paste(parameters[,1])

# filter to just species i need
Species <- c("PBI", "SBI", "SY", "BE", "AH", "POK", "SOK", "ASP", "BPO", "CAR", "ROW", "WWL")

CT <- raster(paste0(dirData,"ESC/ct.tif")) # continentality
DAMS <- raster(paste0(dirData,"ESC/dams.tif")) # detailed aspect method of scoring (wind risk)
SMR <- raster(paste0(dirData,"ESC/smr.tif")) # soil moisture regime
SNR <- raster(paste0(dirData,"ESC/snr.tif")) # soil nutrient regime

# original ESC rasters are at 250m res, CHESS are 1000m - aggregate ESC up to 1000
CT <- aggregate(CT, fact=4,fun=min)
DAMS <- aggregate(DAMS, fact=4,fun=min)
SMR <- aggregate(SMR, fact=4,fun=min)
SNR <- aggregate(SNR, fact=4,fun=min)

stckESC <- stack(CT,DAMS,SMR,SNR)
plot(stckESC)

### function ESC ---------------------------------------------------------------

calculateEscSuitabilityForSpecies <- function(sp,at,ct,da,md,smr,snr,year){
  
  sppar<-as.numeric(subset(parameters,parameters$abbr==sp))
  at4<-sppar[8]; at3<-sppar[9] ;at2<-sppar[10] ;at1<-sppar[11]
  md4<-sppar[12]; md3<-sppar[13] ;md2<-sppar[14] ;md1<-sppar[15]
  ct4<-sppar[16]; ct3<-sppar[17] ;ct2<-sppar[18] ;ct1<-sppar[19]
  da4<-sppar[20]; da3<-sppar[21] ;da2<-sppar[22] ;da1<-sppar[23]
  smr4<-sppar[24]; smr3<-sppar[25] ;smr2<-sppar[26] ;smr1<-sppar[27]
  snr4<-sppar[28]; snr3<-sppar[29] ;snr2<-sppar[30] ;snr1<-sppar[31]
  
  ATF<-at4*at^3 + at3*at^2 + at2*at + at1
  CTF<-ct4*ct^3 + ct3*ct^2 + ct2*ct + ct1
  DAF<-da4*da^3 + da3*da^2 + da2*da + da1
  MDF<-md4*md^3 + md3*md^2 + md2*md + md1
  smrF<-smr4*smr^3 + smr3*smr^2 + smr2*smr + smr1
  snrF<-snr4*snr^3 + snr3*snr^2 + snr2*snr + snr1 
  
  ATF[ATF>1]<-1;   ATF[ATF<0]<-0
  CTF[CTF>1]<-1;   CTF[CTF<0]<-0
  DAF[DAF>1]<-1;   DAF[DAF<0]<-0
  MDF[MDF>1]<-1;   MDF[MDF<0]<-0
  smrF[smrF>1]<-1; smrF[smrF<0]<-0
  snrF[snrF>1]<-1; snrF[snrF<0]<-0
  
  # climate
  # ESC_stack<-stack(CTF,DAF,MDF)
  # SS_min <- stackApply(ESC_stack,indices=c(1,1,1), min, na.rm=F)
  
  # climate and soils
  ESC_stack2<-stack(CTF,DAF,MDF,smrF,snrF)
  SS_min2 <- stackApply(ESC_stack2,indices=c(1,1,1,1,1), min, na.rm=F) 
  
  SUIT<-SS_min2*ATF #remember to change SS_min/2 here depending on if I have soils in or not
  YC<-SUIT*sppar[4]    
  
  projection(SUIT)<-projection(da)
  projection(YC)<-projection(da)
  
  # optional reclassification of suitability for m and m2 yield class, can take it out
  # name1<-paste(sp, names(at),sep = "_")
  # name2<-gsub("AT_","",name1)
  m  <- c(-1, 0.3, 1,  0.3, 0.5, 2,  0.5, 0.75, 3, 0.75, Inf, 4 )
  m2 <- c(-40, 0, 0,  0, 1, 1,  1, 2, 2,  2, 3, 3, 3, 4, 4 , 4, 5, 5, 5, 6, 6,  6, 7, 7,  7, 8, 8, 8, 9, 9 , 9, 10, 10,
          10, 11, 11,  11, 12, 12,  12, 13, 13, 13, 14, 14 , 14, 15, 15, 15, 16, 16,  16, 17, 17,  17, 18, 18, 18, 19, 19 , 19, 20, 20,
          20, 21, 21,  21, 22, 22,  22, 23, 23, 23, 24, 24 , 24, 25, 25, 25, 26, 26,  26, 27, 27,  27, 28, 28, 28, 29, 29 , 29, 30, 30,
          30, 31, 31,  31, 32, 32,  32, 33, 33, 33, 34, 34 , 34, 35, 35, 35, 36, 36,  36, 37, 37,  37, 38, 38, 38, 39, 39 , 39, 40, 40)
  # reclassify data to reduce filesize/simplify visualisation
  suitmat <- matrix(m, ncol=3, byrow=TRUE)
  ycmat   <- matrix(m2, ncol=3, byrow=TRUE)
  
  SUITS <- reclassify(SUIT, suitmat)
  YCS   <- reclassify(YC  , ycmat)
  
  suit  <-paste0( sp,"_soil_suit_",year)
  yc    <-paste0( sp,"_soil_yc_",year)
  
  # don't need suitability for CRAFTY
  #writeRaster(SUITS, filename=paste0("~/Documents/CRAFTY-UK/ESC/",suit), format="GTiff",overwrite=TRUE)
  # baseline path
  #writeRaster(YCS, filename=paste0(dirESC,"/outputs_CHESSbaseline/",yc), format="GTiff",overwrite=TRUE)
  # adjusted CMD
  #writeRaster(YCS, filename=paste0(dirESC,"/outputs_CHESSbaseline/",yc,"_mdAdj"), format="GTiff",overwrite=TRUE)
  # scenarios path
  #writeRaster(YCS, filename=paste0(dirESC,"/outputs_speedRCP85/",yc,"_mdAdj"), format="GTiff",overwrite=TRUE)  
  writeRaster(YCS, filename=paste0(dirESC,"/outputs_speed_future/",yc,"_mdAdj_",rcp), format="GTiff",overwrite=TRUE)
  
}


### Baseline 1991 - 2011 -------------------------------------------------------

AT <- raster(paste0(dirData,"chess_baseline_rst/gdd_1991_2011_annual_rpj.tif"))
MD <- raster(paste0(dirData,"chess_baseline_rst/CMD_adj_1991_2011_baseline_rpj.tif"))

stckESC <- stack(AT,MD,CT,DAMS,SMR,SNR)
plot(stckESC)

for(i in 1:length(Species)){
  calculateEscSuitabilityForSpecies(sp=Species[i],AT, CT, DAMS, MD, SMR, SNR, "baseline")
}

# check
rsts <- list.files(paste0(dirESC,"outputs_CHESSbaseline"),full.names = TRUE)
# non-adjusted
#rsts <- grep("baseline.tif", rsts, value=TRUE)
baselineESC <- do.call(stack, lapply(rsts, raster))
plot(baselineESC)


### Speed future timesteps -----------------------------------------------------

lstRCP <- c("rcp26","rcp45","rcp60","rcp85")

timesteps <- c("2010_2030","2020_2040","2030_2050","2040_2060","2050_2070","2060_2080")

for (rcp in lstRCP){
  
  #rcp <- lstRCP[1]
  
  for (i in timesteps){
    
    #i <- "2010_2030"
    
    AT <- raster(paste0(dirData,"speed_future_rst/gdd_",i,"_annual_",rcp,"_rpj.tif"))
    MD <- raster(paste0(dirData,"speed_future_rst/CMD_adj_",i,"_annual_",rcp,"_rpj.tif"))
    
    for(j in 1:length(Species)){
      calculateEscSuitabilityForSpecies(sp=Species[j],AT, CT, DAMS, MD, SMR, SNR, i)
    }
    
  }
  
}