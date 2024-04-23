
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
lstFolders <- c("chess_baseline","speed_future_rcp26","speed_future_rcp45","speed_future_rcp85")

for (folder in lstFolders){
  
  folder <- lstFolders[1]
  
  lstFiles <- list.files(paste0(dirData,"/",folder))
  
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
    
    rcp <- substring(folder, 14,18)
    
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
  
  write.csv(metadata, file.path(dirData,"/",folder, "metadata.csv"), row.names = FALSE)
  
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
  
  
  ### loop to calculate cmd from pet and pr ----
  
  metadata.filter <- filter(metadata, variable %in% c("pr","precip","pet"))
  
  lstYear <- metadata.filter$from_year
  
  for (folder in lstFolders){
    
    folder <- lstFolders[1]
    
    for (yr in lstYear){
      
      yr <- lstYear[1]
      
      if (folder == "chess_baseline"){
        
        pr <- stars::read_ncdf(paste0(dirData,folder,"/chess-met_precip_gb_1km_20yr-mean-monthly_",yr,"0101-19811231.nc"))
        pr <- st_set_crs(pr, 27700) #  has values
        # convert to monthly total
        pr<- precip_to_total(pr)
        
        pet <- stars::read_ncdf(paste0(dirData,folder,"/chess-pe_pet_uk_1km_20yr-mean-monthly_",yr,"0101-19811231.nc"))
        pet <- st_set_crs(pet, 27700) #  has values
        pet <- pet_to_total(pet)
        
        # calculate monthly moisture deficit (mMD)
        # positive values = deficit, negative values = surplus
        mMD <- pet - pr
        mMD <- st_set_crs(mMD, 27700) # has values
        
        #stars::write_stars(mMD, dsn = paste0(dirScratch,"chess_mMD_1991_2011_monthly.tif"), NA_value = NA)

        # issue here. can't read in tiff as brick, which then affects function below
        #EtoPrDiff <- brick(paste0(dirScratch,"chess_mMD_1991_2011_monthly.tif"))
        
        #EtoPrDiff <- st_apply(mMD[,,3], 1:3, calcCMD, keep = TRUE)
        #plot(EtoPrDiff)
        
        CMD <- st_apply(mMD[,,,1:12], 1:2, calcCMD, .fname = "CMD") # I think this should select times 1:12 (months jan-feb) from mMD, whilst keeping dimension x/y (1:2)
        # I think it bloody works!!!
        plot(CMD)
        
         # apply adjustment
        diff <- 0.0011*CMD^2 - 0.076*CMD + 0.08465
        CMD_adj <- CMD - diff
        
        plot(CMD_adj)
        
        # write as raster, somehow
        write_stars(CMD_adj, paste0(dirScratch,"chess_CMD_annual_",yr,".tif"))

      } else {
        
       # read in speed pattern
        
      }
      
    }
    }
  


### now reproject and extend to same extent as ESC rasters ----

# esc reference raster
reference <- raster(paste0(dirData,"/ESC/ct.tif"))
# aggregate to 1k
reference <- aggregate(reference, fact=4,fun=min)
plot(reference)

files <- list.files(paste0(dirScratch),full.names = T)
files <- Filter(function(x) grepl("CMD|gdd", x), files)
#files <- Filter(function(x) grepl("rcp26|rcp45", x), files)

for (i in files){
  
  i <- files[1]
  
  x <- raster(i)
  # assign projection
  projection(x) <- crs(reference)
  # extend to extent
  x_crop <- extend(x,reference)
  # write reprojected file
  
  file.name <- stringr::str_split(i,pattern = "/")[[1]][5]
  reproj.name <- substr(file.name,1,nchar(file.name)-4)
  
  writeRaster(x_crop, paste0(dirScratch,"/speed_future_rst/",reproj.name,"_rpj.tif"),overwrite=T)
  
  }
  

