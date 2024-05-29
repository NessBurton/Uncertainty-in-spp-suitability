
# date: 28-05-24
# author: VB
# desc: script to assess variability/uncertainty between ESC results for native species

library(terra)

### working dirs ----

wd <- "C:/Users/vbu/OneDrive - the Woodland Trust/Data-analysis/SoWT2_Spp-uncertainty" # WT path
dirData <-  paste0(wd,"/data-raw/")
dirScratch <-  paste0(wd,"/data-scratch/")
dirOut <- paste0(wd,"/data-out/")

### read in ESC results ----

# for each species
lstSpecies <- c("PBI", "SBI", "SY", "BE", "AH", "POK", "SOK", "ASP", "BPO", "CAR", "ROW", "WWL")
#lstScenario <- c("speed_rcp26","speed_rcp45","speed_rcp85")
lstTimesteps <- c("2010-2030","2020-2040","2030-2050","2040-2060","2050-2070","2060-2080")

for (sp in lstSpecies){
  
  #sp <- lstSpecies[12]
  # read in suitability results for all scenarios
  print(paste0("Working on species: ", sp))
  
  # per timestep
  for (ts in lstTimesteps){
    
    #ts <- lstTimesteps[6]
    print(paste0("Timestep: ", ts))
    
    baseline <- rast(paste0(dirOut,"chess_",sp,"_soil_suit_1991-2011.tif"))
    rcp25 <- rast(paste0(dirOut,"speed_rcp26_",sp,"_soil_suit_",ts,".tif"))
    rcp45 <- rast(paste0(dirOut,"speed_rcp45_",sp,"_soil_suit_",ts,".tif"))
    rcp85 <- rast(paste0(dirOut,"speed_rcp85_",sp,"_soil_suit_",ts,".tif"))
    
    # plot them to visually compare first
    par(mfrow= c(2,2)) 
    
    plot(baseline); plot(rcp25); plot(rcp45); plot(rcp85)
    
  }
  
}




### calculate agreement between scenarios a la B4EST work ----