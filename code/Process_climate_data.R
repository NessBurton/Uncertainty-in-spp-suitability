
# date: 16-04-24
# author: VB
# desc: script to process climate data (CHESS baseline 1991-2011 and SPEED future secnarios (downscaled from UKCP18) from CEH), convert nc files to tif, and calculate CMD for ESC

### working dirs ----

#wd <- "~/Documents/Woodland-Trust/Data-Analysis/SoWT2_Spp-uncertainty" # Mac path
wd <- "C:/Users/vbu/OneDrive - the Woodland Trust/Data-analysis/SoWT2_Spp-uncertainty" # WT path
dirData <-  paste0(wd,"/data-raw/")
dirScratch <-  paste0(wd,"/data-scratch/")
dirOut <- paste0(wd,"/data-out/")

### load libs ----

library(dplyr)
library(stars)
library(raster)
library(terra)
library(ncmeta)
library(ncdf4)

# load in functions from #CEH script
source(file.path(wd,"code","NC-to-tiff.R"))

### loop to process all .nc files, per baseline/scenario ----

# list folders within data-raw
lstFolders <- c("chess_baseline","speed_rcp26","speed_rcp45","speed_rcp85")

for (folder in lstFolders){
  
  #folder <- lstFolders[2]
  
  print(paste0("Processing scenario: ",folder))
  
  lstFiles <- list.files(paste0(dirData,folder))
  
  #write.csv(metadata, file.path(dirData,"/",folder, "metadata.csv"), row.names = FALSE)
  
  # process nc files
 
  ### functions ----
  # calculate total monthly potential evapotranspiration
  # the data provided is mean monthly potential evapotranspiration (per day)
  # it needs to be transformed to total - multiply by the number of days in a month
  pet_to_total <- function (x) {
    month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    dims <- dim(x)
    total <- x * array(rep(month_days, each = prod(dims[1:2])), dims)
    total <- st_set_dimensions(total, "time", tolower(month.abb))
    total
  }
  
  # calculate total monthly precipitation
  # the data provided is mean monthly precipitation (per second)
  # needs to be transformed to total - multiply by the number seconds in a month
  precip_to_total <- function (x) {
    month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    month_seconds <- month_days * 3600 * 24
    dims <- dim(pr)
    total <- pr * array(rep(month_seconds, each = prod(dims[1:2])), dims)
    total <- st_set_dimensions(total, "time", tolower(month.abb))
    total
  }
  
  # calculate CMD - function from  CEH
  calcCMD <- function(EtoPrDiff) {
    EtoPrDiff <- EtoPrDiff[!is.na(EtoPrDiff)]
    if (length(EtoPrDiff) > 0) {
      CMD <- 0
      accumulator <- 0
      for (a in EtoPrDiff) {
        if (accumulator >= 0) {
          accumulator <- max(accumulator + a, 0)
          if (accumulator > CMD) CMD <- accumulator
        } else {
          accumulator <- 0
        }
      }
    } else {
      CMD <- NA
    }
    
    return(CMD)
  }
  
  # reproject and extend to same extent as ESC raster
  # esc reference raster
  reference <- rast(paste0(dirData,"/ESC/ct.tif"))
  # aggregate to 1k
  reference <- aggregate(reference, fact=4,fun=min)
  #plot(reference)
  
  
  ### loop to calculate cmd from pet and pr ----
  
  metadata.filter <- filter(metadata, variable %in% c("pr","precip","pet"))
  
  lstFrmYear <- metadata.filter$from_year
  lstToYear <- metadata.filter$to_year
  
  for (folder in lstFolders){
    
    #folder <- lstFolders[2]
    
    if (folder == "chess_baseline"){
      
      # create metadata
      description <- c(
        gdd = "growing degree days and number of days",
        precip = "precipitation",
        pet = "Penman-Monteith potential evapotranspiration, for a well-watered grass surface")
      
      metadata <- tibble::tibble(file = lstFiles) %>%
        tidyr::extract(
          file,
          c("variable", "summary", "temp_resolution", "from", "to"),
          "chess-[[:alpha:]]+_(.+)_[[:alpha:]]{2}_1km_20yr-([[:alpha:]]+)-([[:alpha:]]+)_([[:digit:]]+)-([[:digit:]]+)",
          remove = FALSE) %>%
        mutate(from_year = substr(from, 1, 4),
               from_month = substr(from, 5, 6),
               to_year = substr(to, 1, 4),
               to_month = substr(to, 5, 6),
               description = description[variable])
    } else {
      
      rcp <- substring(folder, 7,11)
      
      # create metadata
      description <- c(
        gdd = "growing degree days and number of days",
        pr = "precipitation",
        pet = "Penman-Monteith potential evapotranspiration, for a well-watered grass surface")
      
      metadata <- tibble::tibble(file = lstFiles) %>%
        tidyr::extract(
          file,
          c("variable", "summary", "temp_resolution", "from", "to"),
          paste0(".+",rcp,"_bias_corrected_01_(.+)_uk_1km_20yr-([[:alpha:]]+)-([[:alpha:]]+)_([[:digit:]]+)-([[:digit:]]+)"),
          remove = FALSE
        ) %>%
        mutate(
          from_year = substr(from, 1, 4),
          from_month = substr(from, 5, 6),
          to_year = substr(to, 1, 4),
          to_month = substr(to, 5, 6),
          description = description[variable]
        )
    }
    
    
    for (i in 1:length(lstFrmYear)){
      
      #i <- 1
      yrFrm <- lstFrmYear[i]
      yrTo <- lstToYear[i]
      
      print(paste0("Processing timestep from: ", yrFrm, " to ", yrTo, " for scenario ", folder))
      
      if (folder == "chess_baseline"){
        
        # read in gdd
        gdd <- stars::read_ncdf(paste0(dirData,folder,"/chess-met_gdd_gb_1km_20yr-mean-annual_",yrFrm,"0101-",yrTo,"1231.nc"))
        gdd <- st_set_crs(gdd, 27700)
        # drop time dimension (annual so only 1)
        gdd <- gdd[,1:656,1:1057,drop = TRUE] #%>%  dim()
        # manually convert to raster
        # reproject and extend to ESC req.
        xyz <- data.frame(x = gdd[,1:656,], y = gdd[,,1:1057], z = gdd["gdd",,])
        gdd2 <- rasterFromXYZ(xyz, crs = 27700)
        gdd2 <- rast(gdd2[[1]], names = c("gdd"))
        #print(plot(gdd2))
        
        # assign projection
        terra::project(gdd2, reference) #) <- crs(reference)
        # extend to extent
        gdd2 <- extend(gdd2,reference)
        #plot(gdd2)
        
        # write to file
        writeRaster(gdd2, paste0(dirScratch,"chess_gdd_",yrFrm,"-",yrTo,"_rpj.tif"),overwrite=T) # seems to be upside down?!
        
        print(paste0("gdd processed & written to file"))
        
        # read in pr
        pr <- stars::read_ncdf(paste0(dirData,folder,"/chess-met_precip_gb_1km_20yr-mean-monthly_",yrFrm,"0101-",yrTo,"1231.nc"))
        pr <- st_set_crs(pr, 27700) #  has values
        # convert to monthly total
        pr<- precip_to_total(pr)
        
        pet <- stars::read_ncdf(paste0(dirData,folder,"/chess-pe_pet_uk_1km_20yr-mean-monthly_",yrFrm,"0101-",yrTo,"1231.nc"))
        pet <- st_set_crs(pet, 27700) #  has values
        pet <- pet_to_total(pet)
        
        # calculate monthly moisture deficit (mMD)
        # positive values = deficit, negative values = surplus
        mMD <- pet - pr
        mMD <- st_set_crs(mMD, 27700) # has values
        
        #calculate climatic moisture deficit (CMD)
        CMD <- st_apply(mMD[,,,1:12], 1:2, calcCMD, .fname = "CMD") # I think this should select times 1:12 (months jan-feb) from mMD, whilst keeping dimension x/y (1:2)
        plot(CMD)
        
        # apply adjustment ( to account for penman montieth version of pet)
        diff <- 0.0011*CMD^2 - 0.076*CMD + 0.08465
        CMD_adj <- CMD - diff
        
        print(plot(CMD_adj))
        
        # manually extract values to get to a raster
        xyz <- data.frame(x = CMD_adj[,1:656,], y = CMD_adj[,,1:1057], z = CMD_adj["CMD",,])
        CMD2 <- rasterFromXYZ(xyz, crs = 27700)
        CMD2 <- rast(CMD2[[1]], names = c("CMD"))
        print(plot(CMD2))
   
        # assign projection
        terra::project(CMD2, reference) #) <- crs(reference)
        # extend to extent
        CMD2 <- extend(CMD2,reference)
        
        writeRaster(CMD2, paste0(dirScratch,"chess_CMD_",yrFrm,"-",yrTo,"_rpj.tif"), overwrite=T)
        
        print(paste0("CMD calculated and written to file"))

      } else {
        
        # read in speed pattern
        
        # read in gdd
        gdd <- stars::read_ncdf(paste0(dirData,folder,"/ukcp18-",folder,"_bias_corrected_01_gdd_uk_1km_20yr-mean-annual_",yrFrm,"1201-",yrTo,"1130.nc"))
        gdd <- st_set_crs(gdd, 27700)
        # drop time dimension (annual so only 1)
        gdd <- gdd[,1:656,1:1057,drop = TRUE] #%>%  dim()
        # manually convert to raster
        xyz <- data.frame(x = gdd[,1:656,], y = gdd[,,1:1057], z = gdd["gdd",,])
        gdd2 <- rasterFromXYZ(xyz, crs = 27700)
        gdd2 <- rast(gdd2[[1]], names = c("gdd"))
        print(plot(gdd2))
        
        # reproject and extend to ESC req.
        # assign projection
        terra::project(gdd2, reference) #) <- crs(reference)
        # extend to extent
        gdd2 <- extend(gdd2,reference)
        #plot(gdd2)
        
        # write to file
        writeRaster(gdd2, paste0(dirScratch,"speed_gdd_",yrFrm,"-",yrTo,"_rpj.tif"),overwrite=T) # seems to be upside down?!
        
        print(paste0("gdd processed & written to file"))
        
        # read in pr
        pr <- stars::read_ncdf(paste0(dirData,folder,"/ukcp18-",folder,"_bias_corrected_01_pr_uk_1km_20yr-mean-monthly_",yrFrm,"12-",yrTo,"11.nc"))
        pr <- st_set_crs(pr, 27700) #  has values
        # convert to monthly total
        pr<- precip_to_total(pr)
        
        # read in pet
        pet <- stars::read_ncdf(paste0(dirData,folder,"/ukcp18-",folder,"_bias_corrected_01_pet_uk_1km_20yr-mean-monthly_",yrFrm,"12-",yrTo,"11.nc"))
        pet <- st_set_crs(pet, 27700) #  has values
        pet <- pet_to_total(pet)
        
        # calculate monthly moisture deficit (mMD)
        # positive values = deficit, negative values = surplus
        mMD <- pet - pr
        mMD <- st_set_crs(mMD, 27700) # has values
        
        # calculate climatic moisture deficit (CMD)
        CMD <- st_apply(mMD[,,,1:12], 1:2, calcCMD, .fname = "CMD") # I think this should select times 1:12 (months jan-feb) from mMD, whilst keeping dimension x/y (1:2)
        plot(CMD)
        
        # apply adjustment ( to account for penman montieth version of pet)
        diff <- 0.0011*CMD^2 - 0.076*CMD + 0.08465
        CMD_adj <- CMD - diff
        
        plot(CMD_adj)
        
        # manually extract values to get to a raster
        xyz <- data.frame(x = CMD_adj[,1:656,], y = CMD_adj[,,1:1057], z = CMD_adj["CMD",,])
        CMD2 <- rasterFromXYZ(xyz, crs = 27700)
        CMD2 <- rast(CMD2[[1]], names = c("CMD"))
        #print(plot(CMD2))
        
        # assign projection
        terra::project(CMD2, reference) #) <- crs(reference)
        # extend to extent
        CMD2 <- extend(CMD2,reference)
        plot(CMD2)
        
        writeRaster(CMD2, paste0(dirScratch,"speed_CMD_",yrFrm,"-",yrTo,"_rpj.tif"), overwrite=T)
        
        print(paste0("CMD calculated and written to file"))
        
          
        }
        
      }
      
    }
  }

