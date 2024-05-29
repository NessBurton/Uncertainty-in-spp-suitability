
# date: 28-05-24
# author: VB
# desc: script to assess variability/uncertainty between ESC results for native species

library(terra)
library(tidyverse)
library(ggplot2)
library(viridis)

### working dirs ----

wd <- "C:/Users/vbu/OneDrive - the Woodland Trust/Data-analysis/SoWT2_Spp-uncertainty" # WT path
dirData <-  paste0(wd,"/data-raw/")
dirScratch <-  paste0(wd,"/data-scratch/")
dirOut <- paste0(wd,"/data-out/")
dirFigs <- paste0(wd,"/figures/")


### read in ESC results ----

# for each species
lstSpecies <- c("PBI", "SBI", "SY", "BE", "AH", "POK", "SOK", "ASP", "BPO", "CAR", "ROW", "WWL")
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

# list all tifs
rsts <-  list.files(paste0(dirOut),pattern = "*.tif", full.names = T)

# list all timesteps
lstTimesteps <- c("2010-2030","2020-2040","2030-2050","2040-2060","2050-2070","2060-2080")
# list all species
lstSpecies <- c("PBI", "SBI", "SY", "BE", "AH", "POK", "SOK", "ASP", "BPO", "CAR", "ROW", "WWL")

# for each timestep
for (ts in lstTimesteps){
  
  #ts <- lstTimesteps[2]
  # filter to timestep
  rstsTS <- grep(ts, rsts, value = TRUE)
  # filter to just suitability, not yield
  rstsTS <- grep(c("speed"), rstsTS, value = T)
  rstsTS <- grep(c("suit"), rstsTS, value = T)
  
  for (sp in lstSpecies){
    
    #sp <- lstSpecies[1]
    
    SP.name <- ifelse(grepl("PBI", sp), 'Downy birch', 
                      ifelse(grepl("SBI", sp), 'Silver birch',
                             ifelse(grepl("SY", sp), 'Sycamore',
                                    ifelse(grepl("BE", sp), 'Beech',
                                           ifelse(grepl("AH", sp), 'Ash',
                                                  ifelse(grepl("POK", sp), 'Pedunculate oak',
                                                         ifelse(grepl("SOK", sp), 'Sessile oak',
                                                                ifelse(grepl("ASP", sp), 'Aspen',
                                                                       ifelse(grepl("BPO", sp), "Black poplar",
                                                                              ifelse(grepl("CAR", sp), "Alder",
                                                                                     ifelse(grepl("ROW", sp), "Rowan",
                                                                                            ifelse(grepl("WWL", sp), "White willow", NA))))))))))))
    
    # filter to species
    rstsSpp <- grep(sp, rstsTS, value = TRUE)
    
    # raster stack
    rstStack <- rast(rstsSpp)
    #plot(rstStack)
    
    # reclass
    r <- c(0, 1, NA,
           1, 2, 1,
           2, 3, 1,
           3, 4, 1,
           4, 5, 1)
    
    r_matrix <- matrix(r, ncol=3, byrow=TRUE)
    
    rstReclass <- classify(rstStack, r_matrix)
    plot(rstReclass)
    
    # sum to find agreement
    sumStack <- sum(rstReclass, na.rm = T)
    plot(sumStack)
    
    # Convert raster to dataframe
    df1 <- as.data.frame(sumStack, xy=T)
    names(df1) <- c("x", "y", "RCPagree")
    
    (p1 <- ggplot(data = df1) + 
        geom_tile(data = df1 %>% filter(!is.na(RCPagree)), mapping = aes(x = x, y = y, fill = RCPagree), size = 1)+
        scale_fill_viridis(direction = -1)+
        theme_bw()+
        ggtitle(paste0(SP.name," | ",ts))+
        theme(plot.title = element_text(face="bold",size=22),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()))
    
    ggsave(p1, file=paste0(dirFigs,"RCP_agreement_",sp,"_timestep_",ts,".png"), width=8, height=10, dpi=300)
    
  }
  
}
