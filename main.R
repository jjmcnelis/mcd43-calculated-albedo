#####################################################################################
# Set up directories
#####################################################################################

wkd <- 'D:/Projects/MCD43A/github/'
setwd(wkd)

data_dir <- paste(wkd, 'example-data/', sep='') #HOME
out_dir <- paste(data_dir, 'output/', sep = '')
suppressWarnings(dir.create(out_dir))

# Solar Optical Depth (Value between 0.00 and 0.98 x0.02 with 2-digit precision)
SolOptDepth <- 0.20

#####################################################################################
# Import libraries
#####################################################################################

library(ncdf4)
library(data.table)
library(RAtmosphere)
library(raster)

#####################################################################################
# Open helper files
#####################################################################################

# Open look up tables
vislut <- t(read.csv('vis_lut.csv'))
nirlut <- t(read.csv('nir_lut.csv'))
swlut <- t(read.csv('sw_lut.csv'))
# External functions
source('f_albedo.R')
source('proc_albedo.R')

#####################################################################################
# NetCDF ops
#####################################################################################

nc_in = paste(data_dir, 'MCD43A1.006_aid0001_large.nc', sep = '')
nc <- nc_open(nc_in)

# Print variables
print(attributes(nc$var)$names)
# Print dimensions
print(attributes(nc$dim)$names)

# Get data.frame of latlon combinations and separate to vectors
nc_latlon <- expand.grid(nc$dim$lat$vals, nc$dim$lon$vals)
lat <- nc_latlon[1]
lon <- nc_latlon[2]
# Get CRS
crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
# Get time dimension values as array and length of time dimension
nc_time <- nc$dim$time$vals
nt <- length(nc_time)

print(paste('There are',nrow(nc_latlon), 'pixels and', 
            length(nc_time), 'timesteps in this dataset.', sep = ' '))

#####################################################################################
# Process visible, near-infrared, and shortwave albedo (proc_albedo.R, f_albedo.R)
# General steps:
# Mask QC
# Calculate black sky albedo
# Calculate white sky albedo
# Calculate actual albedo
#####################################################################################

vis_data <- proc_albedo('vis')
nir_data <- proc_albedo('nir')
sw_data <- proc_albedo('shortwave')

#####################################################################################
# Write new NetCDF file
#####################################################################################

# Change RasterStacks to matrices because ncdf4 package doesn't seem to like them
vis_black <- as.matrix(vis_data[[1]])
vis_white <- as.matrix(vis_data[[2]])
vis_actual <- as.matrix(vis_data[[3]])
nir_black <- as.matrix(nir_data[[1]])
nir_white <- as.matrix(nir_data[[2]])
nir_actual <- as.matrix(nir_data[[3]])
shortwave_black <- as.matrix(sw_data[[1]])
shortwave_white <- as.matrix(sw_data[[2]])
shortwave_actual <- as.matrix(sw_data[[3]])

# Define output dimensions
outlat <- ncdim_def('lat', nc$dim$lat$units, as.double(nc$dim$lat$vals), create_dimvar=TRUE)
outlon <- ncdim_def('lon', nc$dim$lon$units, as.double(nc$dim$lon$vals), create_dimvar=TRUE)
outtime <- ncdim_def('time', nc$dim$time$units, as.integer(nc$dim$time$vals), create_dimvar=TRUE, unlim = TRUE)

# Define new variables for calculated albedo
vis_bsa <- ncvar_def('visible_black_sky_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
vis_wsa <- ncvar_def('visible_white_sky_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
vis_asa <- ncvar_def('visible_actual_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
nir_bsa <- ncvar_def('nir_black_sky_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
nir_wsa <- ncvar_def('nir_white_sky_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
nir_asa <- ncvar_def('nir_actual_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
sw_bsa <- ncvar_def('shortwave_black_sky_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
sw_wsa <- ncvar_def('shortwave_white_sky_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")
sw_asa <- ncvar_def('shortwave_actual_albedo', units = '', dim = list(outlat, outlon, outtime), missval = -9999, prec = "float")

# Set output filename and create output netcdf file
albedo_out <- paste0(out_dir, 'MCD43A_calcalbedo.nc')
nc_out <- nc_create(albedo_out, list(vis_bsa, vis_wsa, vis_asa, nir_bsa, nir_wsa, nir_asa, sw_bsa, sw_wsa, sw_asa))

# Add data to new variables
for (tstep in 1:nt){
  print(paste("Writing timestep", tstep, "to netcdf variables . . .", sep = " "))
  ncvar_put(nc_out, vis_bsa, vis_black[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, vis_wsa, vis_white[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, vis_asa, vis_actual[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, nir_bsa, nir_black[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, nir_wsa, nir_white[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, nir_asa, nir_actual[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, sw_bsa, shortwave_black[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, sw_wsa, shortwave_white[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, sw_asa, shortwave_actual[,tstep], start = c(1,1,tstep), count = c(-1, -1, 1))
}

# Write out and close
nc_out <- nc_open(albedo_out, write = TRUE)
nc_close(nc_out)
