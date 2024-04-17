
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
  
  folder <- lstFolders[4]
  
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
  
  # loop to read each file, convert from averages to totals if either pet or pr, and write to tif
  
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
    dims <- dim(x)
    total <- x * array(rep(month_seconds, each = prod(dims[1:2])), dims)
    total <- st_set_dimensions(total, "time", tolower(month.abb))
    total
  }
  
  
  for (i in 1:nrow(metadata)){
    
    #i <- 1
    
    #x <- read_nc(paste0(dirData, metadata$file[i]))
    x <- stars::read_ncdf(file.path(dirData,folder,"/",metadata$file[i]), 27700)
    
    # if pet file, convert monthly to total
    if (metadata$variable[i] == "pet"){
      
      x <- pet_to_total(x)
      
    }
    
    # if pr file, convert monthly to total
    if (metadata$variable[i] %in% c("pr","precip")){
      
      x <- precip_to_total(x)
      
    }
    
    if (folder == "chess_baseline"){
      
      stars::write_stars(x, paste0(dirScratch,"chess_",metadata$variable[i],"_",metadata$from_year[i],"_",metadata$to_year[i],"_",metadata$temp_resolution[i],".tif"))
      
    } else {
      
      stars::write_stars(x, dsn = paste0(dirScratch,"speed_",metadata$variable[i],"_",metadata$from_year[i],"_",metadata$to_year[i],"_",metadata$temp_resolution[i],"_",rcp,".tif"))

    }
    
  }
  
  # now calculate CMD = climatic moisture deficit, required for ESC
  # use monthly values for pet and prec
  
  # function
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
  
  # need to do this per timestep
  timesteps <- c("2010_2030","2020_2040","2030_2050","2040_2060","2050_2070","2060_2080")
  
  for (i in timesteps){
    
    if (folder == "chess_baseline"){
      
      i <- "1991_2011"
      
      pet <- read_stars(paste0(dirScratch,"/chess_pet_",i,"_monthly.tif"))
      prec <- read_stars(paste0(dirScratch,"/chess_precip_",i,"_monthly.tif"))
      
    } else {
      
      pet <- read_stars(paste0(dirScratch,"/speed_pet_",i,"_monthly_",rcp,".tif"))
      prec <- read_stars(paste0(dirScratch,"/speed_pr_",i,"_monthly_",rcp,".tif"))
      
      }
    
    # calculate monthly moisture deficit (mMD)
    # positive values = deficit, negative values = surplus
    mMD <- pet - prec
  
    #dev.off()
    #plot(mMD, key.pos = 4, key.width = lcm(1), box_col = "white", col = hcl.colors(10))
    
    #stars::write_stars(mMD, dsn = paste0(dirScratch,"chess_mMD_1991_2011_monthly.tif"), NA_value = NA)
    
    # issue here. can't read in tiff as brick, which then affects function below
    #EtoPrDiff <- brick(paste0(dirScratch,"chess_mMD_1991_2011_monthly.tif"))
    
    
    # work around - subset each month from stars object and convert to raster
    jan <- raster(mMD[[1]][,,1], crs = 27700)
    feb <- raster(mMD[[1]][,,2], crs = 27700)
    mar <- raster(mMD[[1]][,,3], crs = 27700)
    apr <- raster(mMD[[1]][,,4], crs = 27700)
    may <- raster(mMD[[1]][,,5], crs = 27700)
    jun <- raster(mMD[[1]][,,6], crs = 27700)
    jul <- raster(mMD[[1]][,,7], crs = 27700)
    aug <- raster(mMD[[1]][,,8], crs = 27700)
    sep <- raster(mMD[[1]][,,9], crs = 27700)
    oct <- raster(mMD[[1]][,,10], crs = 27700)
    nov <- raster(mMD[[1]][,,11], crs = 27700)
    dec <- raster(mMD[[1]][,,12], crs = 27700)
    
    EtoPrDiff <- brick(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,dec)
  
    rasterCMD <- raster::calc(EtoPrDiff, fun = calcCMD)
    #rasterCMD <- terra::app(EtoPrDiff, fun = calcCMD)
  
    #dev.off()
    #plot(rasterCMD, col = hcl.colors(10))
  
    # apply adjustment
    diff <- 0.0011*rasterCMD^2 - 0.076*rasterCMD + 0.08465
  
    rasterCMD_adj <- rasterCMD - diff
  
    #dev.off()
    #plot(rasterCMD_adj, col = hcl.colors(10))
  
    #writeRaster(rasterCMD_adj, paste0(dirScratch,"/chess_CMD_adj_1991_2011_baseline.tif") , overwrite=TRUE)
  
    if (folder == "chess_baseline"){
    
      writeRaster(rasterCMD_adj, paste0(dirScratch,"chess_CMD_adj_1991_2011.tif"), overwrite = TRUE)
    
      } else {
    
        writeRaster(rasterCMD_adj, paste0(dirScratch,"speed",rcp,"_CMD_adj_",i,".tif"), overwrite = TRUE)
    
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
  

