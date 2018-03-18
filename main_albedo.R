#####################################################################################
# Set up directories and solar optical depth << change these >>
#####################################################################################

# Set your working directory. Contains R scripts (main_aledo.R, f_albedo.R, proc_albedo.R), lookup tables (vis_lut.csv, nir_lut, sw_lut.csv), and input NetCDF
wkd <- '/path/to/your/working/directory/'
setwd(wkd)

# Set path and filename for input NetCDF. Default output from AppEEARS should be MCD43A1.006_aid0001.nc
nc_in = paste0(wkd, 'example-data/', 'MCD43A1.006_aid0001.nc')

# Set solar optical depth (value between 0.00 and 0.98 x0.02 with 2-digit precision)
SolOptDepth <- 0.20
# Set solar zenith angle (0.0-89.0 degrees, or "local" for local solar noon)
SolZenAngle <- 'local'

# Set path to output NetCDF
out_dir <- paste0(wkd, 'example-output/')
out_filename <- 'MCD43A_calcalbedo.nc'
suppressWarnings(dir.create(out_dir))

#####################################################################################
# Import libraries
#####################################################################################

library(ncdf4)
library(raster)
library(RAtmosphere)

#####################################################################################
# Read external files and functions
#####################################################################################

# Open look up tables
vislut <- t(read.csv('vis_lut.csv'))
nirlut <- t(read.csv('nir_lut.csv'))
swlut <- t(read.csv('sw_lut.csv'))
# External functions
source('f_albedo.R')

#####################################################################################
# Read NetCDF; get variables and dimensions
#####################################################################################

nc <- nc_open(nc_in)

# Get data.frame of latlon combinations and separate to vectors
nc_latlon <- expand.grid(nc$dim$lat$vals, nc$dim$lon$vals)
lat <- nc_latlon[1]
lon <- nc_latlon[2]

# Set coordinate reference system
crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

# Get time dimension values as array and number of timesteps
nc_time <- nc$dim$time$vals
nt <- length(nc_time)

print(paste('There are',nrow(nc_latlon), 'pixels and', 
            length(nc_time), 'timesteps in this dataset.', sep = ' '))

# Define output dimensions
outlat <- ncdim_def('lat', nc$dim$lat$units, as.double(nc$dim$lat$vals), create_dimvar=TRUE)
outlon <- ncdim_def('lon', nc$dim$lon$units, as.double(nc$dim$lon$vals), create_dimvar=TRUE)
outtime <- ncdim_def('time', nc$dim$time$units, as.integer(nc_time), create_dimvar=TRUE, unlim = TRUE)

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

# Get CRS variable from input NetCDF and write put attributes
crs_var <- nc$var[['crs']]

# Set output filename and create output netcdf file
albedo_out <- paste0(out_dir, out_filename)
nc_out <- nc_create(albedo_out, list(vis_bsa, vis_wsa, vis_asa, nir_bsa, nir_wsa, nir_asa, sw_bsa, sw_wsa, sw_asa, crs_var))

#####################################################################################
# Process visible, near-infrared, and shortwave albedo (f_albedo.R)
# 1. Mask QC
# 2. Calculate black sky albedo
# 3. Calculate white sky albedo
# 4. Calculate actual albedo
#####################################################################################

###################################### Visible ######################################

# Get variable
var <- nc$var[['BRDF_Albedo_Parameters_vis']]
# Get QC variable
varqc <- nc$var[['BRDF_Albedo_Band_Mandatory_Quality_vis']]

# Get _FillValue
fillvalue <- ncatt_get(nc, 'BRDF_Albedo_Parameters_vis', "_FillValue")

# Get variable size and dimension count (4)
varsize <- var$varsize
varndims <- var$ndims

# Get qc size and dimension count (3)
qcndims <- varqc$ndims
qcsize <- varqc$varsize

# Loop through timesteps
for(tstep in 1:nt){
  
  print(paste("Processing visible band timestep", tstep, "of", nt, ". . .", sep = " "))
  
  # Convert tstep to real date
  dt <- as.Date(nc_time[tstep], origin = as.Date("2000-01-01"))
  
  # Calculate array of solar zenith numbers for each pixel for timestep n
  if(SolZenAngle == 'local'){
    SZN <- (mapply(SZA, timein = rep(dt, nrow(nc_latlon)), Lat = nc_latlon[,1], Lon = nc_latlon[,2]) - 90)
  }
  if(SolZenAngle != 'local'){
    SZN <- rep(SolZenAngle, nrow(nc_latlon))
  }
  
  # Get array of lookup table values for each solar zenith number
  LUTarr <- as.vector(mapply(LUT, band = rep('vis', length(SZN)), szn = SZN, sod = rep(SolOptDepth, length(SZN))))
  
  # Initialize start and count to read one timestep of the variable
  vstart <- rep(1,varndims)
  vstart[varndims] <- tstep
  vcount <- varsize
  vcount[varndims] <- 1
  
  # Initialize start and count to read one timestep of the variable qc data
  qstart <- rep(1,qcndims)
  qstart[qcndims] <- tstep
  qcount <- qcsize
  qcount[qcndims] <- 1
  
  # Read variable data for timestep n from nc structure
  vardata <- ncvar_get(nc, var, start=vstart, count=vcount)
  # Read qc data for timestep n from nc structure
  qcdata <- ncvar_get(nc, varqc, start=qstart, count=qcount)
  
  # Set _FillValue to NA
  vardata[vardata == fillvalue$value] <- NA
  
  # Get 3 parameters as raster objects
  vardata1 <- raster(t(vardata[1,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  vardata2 <- raster(t(vardata[2,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  vardata3 <- raster(t(vardata[3,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  
  # Apply QC mask to vardata (valid == 0 or 1)
  vardata[!qcdata %in% c(0,1)] <- NA
  
  # Calculate black sky albedo for timestep n
  nBlack <- BlackSA(vardata1, vardata2, vardata3, SZN)
  
  # Calculate white sky albedo for timestep n
  nWhite <- WhiteSA(vardata1, vardata2, vardata3)
  
  # Calculate actual albedo for timestep n
  nActual <- ActualSA(nWhite, nBlack, LUTarr)
  
  ncvar_put(nc_out, vis_bsa, as.matrix(nBlack), start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, vis_wsa, as.matrix(nWhite), start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, vis_asa, as.matrix(nActual), start = c(1,1,tstep), count = c(-1, -1, 1))
  
}

# Clear temporary variables from workspace
remove(var, varqc, fillvalue, varsize, varndims, qcndims, qcsize, dt, SZN, LUTarr, vstart, vcount,
       qstart, qcount, vardata, qcdata, vardata1, vardata2, vardata3, nBlack, nWhite, nActual)

################################### Near-infrared ###################################

# Get variable
var <- nc$var[['BRDF_Albedo_Parameters_nir']]
# Get QC variable
varqc <- nc$var[['BRDF_Albedo_Band_Mandatory_Quality_nir']]

# Get _FillValue
fillvalue <- ncatt_get(nc, 'BRDF_Albedo_Parameters_nir', "_FillValue")

# Get variable size and dimension count (4)
varsize <- var$varsize
varndims <- var$ndims

# Get qc size and dimension count (3)
qcndims <- varqc$ndims
qcsize <- varqc$varsize

# Loop through timesteps
for(tstep in 1:nt){
  
  print(paste("Processing near-infrared band timestep", tstep, "of", nt, ". . .", sep = " "))
  
  # Convert tstep to real date
  dt <- as.Date(nc_time[tstep], origin = as.Date("2000-01-01"))
  
  # Calculate array of solar zenith numbers for each pixel for timestep n
  #################################################################################
  # ORNL DAAC uses an executable (32-bit) that takes date and latitude as input and
  # returns a "solar zenith number" between 0-90. Is this the same as SZA?
  # Need to figure out how that calculation works. For now, results will not be the
  # same as MODIS Global Tool output.
  #################################################################################
  if(SolZenAngle == 'local'){
    SZN <- (mapply(SZA, timein = rep(dt, nrow(nc_latlon)), Lat = nc_latlon[,1], Lon = nc_latlon[,2]) - 90)
  }
  if(SolZenAngle != 'local'){
    SZN <- rep(SolZenAngle, nrow(nc_latlon))
  }
  
  # Get array of lookup table values for each solar zenith number
  LUTarr <- as.vector(mapply(LUT, band = rep('nir', length(SZN)), szn = SZN, sod = rep(SolOptDepth, length(SZN))))
  
  # Initialize start and count to read one timestep of the variable
  vstart <- rep(1,varndims)
  vstart[varndims] <- tstep
  vcount <- varsize
  vcount[varndims] <- 1
  
  # Initialize start and count to read one timestep of the variable qc data
  qstart <- rep(1,qcndims)
  qstart[qcndims] <- tstep
  qcount <- qcsize
  qcount[qcndims] <- 1
  
  # Read variable data for timestep n from nc structure
  vardata <- ncvar_get(nc, var, start=vstart, count=vcount)
  # Read qc data for timestep n from nc structure
  qcdata <- ncvar_get(nc, varqc, start=qstart, count=qcount)
  
  # Set _FillValue to NA
  vardata[vardata == fillvalue$value] <- NA
  
  # Get 3 parameters as raster objects
  vardata1 <- raster(t(vardata[1,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  vardata2 <- raster(t(vardata[2,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  vardata3 <- raster(t(vardata[3,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  
  # Apply QC mask to vardata (valid == 0 or 1)
  vardata[!qcdata %in% c(0,1)] <- NA
  
  # Calculate black sky albedo for timestep n
  nBlack <- BlackSA(vardata1, vardata2, vardata3, SZN)
  
  # Calculate white sky albedo for timestep n
  nWhite <- WhiteSA(vardata1, vardata2, vardata3)
  
  # Calculate actual albedo for timestep n
  nActual <- ActualSA(nWhite, nBlack, LUTarr)
  
  ncvar_put(nc_out, nir_bsa, as.matrix(nBlack), start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, nir_wsa, as.matrix(nWhite), start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, nir_asa, as.matrix(nActual), start = c(1,1,tstep), count = c(-1, -1, 1))
  
}

# Clear temporary variables from workspace
remove(var, varqc, fillvalue, varsize, varndims, qcndims, qcsize, dt, SZN, LUTarr, vstart, vcount,
       qstart, qcount, vardata, qcdata, vardata1, vardata2, vardata3, nBlack, nWhite, nActual)

##################################### Shortwave #####################################

# Get variable
var <- nc$var[['BRDF_Albedo_Parameters_shortwave']]
# Get QC variable
varqc <- nc$var[['BRDF_Albedo_Band_Mandatory_Quality_shortwave']]

# Get _FillValue
fillvalue <- ncatt_get(nc, 'BRDF_Albedo_Parameters_shortwave', "_FillValue")

# Get variable size and dimension count (4)
varsize <- var$varsize
varndims <- var$ndims

# Get qc size and dimension count (3)
qcndims <- varqc$ndims
qcsize <- varqc$varsize

# Loop through timesteps
for(tstep in 1:nt){
  
  print(paste("Processing shortwave infrared band timestep", tstep, "of", nt, ". . .", sep = " "))
  
  # Convert tstep to real date
  dt <- as.Date(nc_time[tstep], origin = as.Date("2000-01-01"))
  
  # Calculate array of solar zenith numbers for each pixel for timestep n
  #################################################################################
  # ORNL DAAC uses an executable (32-bit) that takes date and latitude as input and
  # returns a "solar zenith number" between 0-90. Is this the same as SZA?
  # Need to figure out how that calculation works. For now, results will not be the
  # same as MODIS Global Tool output.
  #################################################################################
  if(SolZenAngle == 'local'){
    SZN <- (mapply(SZA, timein = rep(dt, nrow(nc_latlon)), Lat = nc_latlon[,1], Lon = nc_latlon[,2]) - 90)
  }
  if(SolZenAngle != 'local'){
    SZN <- rep(SolZenAngle, nrow(nc_latlon))
  }
  
  # Get array of lookup table values for each solar zenith number
  LUTarr <- as.vector(mapply(LUT, band = rep('shortwave', length(SZN)), szn = SZN, sod = rep(SolOptDepth, length(SZN))))
  
  # Initialize start and count to read one timestep of the variable
  vstart <- rep(1,varndims)
  vstart[varndims] <- tstep
  vcount <- varsize
  vcount[varndims] <- 1
  
  # Initialize start and count to read one timestep of the variable qc data
  qstart <- rep(1,qcndims)
  qstart[qcndims] <- tstep
  qcount <- qcsize
  qcount[qcndims] <- 1
  
  # Read variable data for timestep n from nc structure
  vardata <- ncvar_get(nc, var, start=vstart, count=vcount)
  # Read qc data for timestep n from nc structure
  qcdata <- ncvar_get(nc, varqc, start=qstart, count=qcount)
  
  # Set _FillValue to NA
  vardata[vardata == fillvalue$value] <- NA
  
  # Get 3 parameters as raster objects
  vardata1 <- raster(t(vardata[1,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  vardata2 <- raster(t(vardata[2,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  vardata3 <- raster(t(vardata[3,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  
  # Apply QC mask to vardata (valid == 0 or 1)
  vardata[!qcdata %in% c(0,1)] <- NA
  
  # Calculate black sky albedo for timestep n
  nBlack <- BlackSA(vardata1, vardata2, vardata3, SZN)
  
  # Calculate white sky albedo for timestep n
  nWhite <- WhiteSA(vardata1, vardata2, vardata3)
  
  # Calculate actual albedo for timestep n
  nActual <- ActualSA(nWhite, nBlack, LUTarr)
  
  ncvar_put(nc_out, sw_bsa, as.matrix(nBlack), start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, sw_wsa, as.matrix(nWhite), start = c(1,1,tstep), count = c(-1, -1, 1))
  ncvar_put(nc_out, sw_asa, as.matrix(nActual), start = c(1,1,tstep), count = c(-1, -1, 1))
  
}

# Clear temporary variables from workspace
remove(var, varqc, fillvalue, varsize, varndims, qcndims, qcsize, dt, SZN, LUTarr, vstart, vcount,
       qstart, qcount, vardata, qcdata, vardata1, vardata2, vardata3, nBlack, nWhite, nActual)

#####################################################################################
# Set attributes in NetCDF variables
#####################################################################################

# Put attributes in crs variable
ncatt_put(nc_out, 'crs', attname = 'epsg_code', attval = 4326, prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'crs', attname = 'horizontal_datum_name', attval = 'WGS84', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'crs', attname = 'semi_major_axis', attval = 6378137, prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'crs', attname = 'inverse_flattening', attval = 298.257223563, prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'crs', attname = 'longitude_of_prime_meridian', attval = 0.0, prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'crs', attname = 'grid_mapping_name', attval = 'latitude_longitude', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'crs', attname = '_CoordinateAxisTypes', attval = 'GeoX GeoY', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'crs', attname = '_CoordinateTransformType', attval = 'Projection', prec = NA, verbose=FALSE, definemode=FALSE )

# Put attributes in visible_*_albedo variables
ncatt_put(nc_out, 'visible_black_sky_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_white_sky_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_actual_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_black_sky_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_white_sky_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_actual_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_black_sky_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_white_sky_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_actual_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_black_sky_albedo', attname = 'long_name', attval = 'Visible_black_sky_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_white_sky_albedo', attname = 'long_name', attval = 'Visible_white_sky_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_actual_albedo', attname = 'long_name', attval = 'Visible_actual_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_black_sky_albedo', attname = 'Description', attval = 'Visible band (8) black sky albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_white_sky_albedo', attname = 'Description', attval = 'Visible band (8) white sky albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'visible_actual_albedo', attname = 'Description', attval = 'Visible band (8) actual (blue sky) albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )

# Put attributes in nir_*_albedo variables
ncatt_put(nc_out, 'nir_black_sky_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_white_sky_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_actual_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_black_sky_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_white_sky_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_actual_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_black_sky_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_white_sky_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_actual_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_black_sky_albedo', attname = 'long_name', attval = 'NIR_black_sky_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_white_sky_albedo', attname = 'long_name', attval = 'NIR_white_sky_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_actual_albedo', attname = 'long_name', attval = 'NIR_actual_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_black_sky_albedo', attname = 'Description', attval = 'Near-infrared band (9) black sky albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_white_sky_albedo', attname = 'Description', attval = 'Near-infrared band (9) white sky albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'nir_actual_albedo', attname = 'Description', attval = 'Near-infrared band (9) actual (blue sky) albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )

# Put attributes in shortwave_*_albedo variables
ncatt_put(nc_out, 'shortwave_black_sky_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_white_sky_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_actual_albedo', attname = 'coordinates', attval = 'time lat lon', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_black_sky_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_white_sky_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_actual_albedo', attname = 'grid_mapping', attval = 'crs', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_black_sky_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_white_sky_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_actual_albedo', attname = 'units', attval = 'unitless', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_black_sky_albedo', attname = 'long_name', attval = 'Shortwave_black_sky_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_white_sky_albedo', attname = 'long_name', attval = 'Shortwave_white_sky_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_actual_albedo', attname = 'long_name', attval = 'Shortwave_actual_albedo', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_black_sky_albedo', attname = 'Description', attval = 'Shortwave-infrared band (10) black sky albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_white_sky_albedo', attname = 'Description', attval = 'Shortwave-infrared band (10) white sky albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'shortwave_actual_albedo', attname = 'Description', attval = 'Shortwave-infrared band (10) actual (blue sky) albedo calculated from MCD43A1', prec = NA, verbose=FALSE, definemode=FALSE )

# Put attributes in lat and lon variables
ncatt_put(nc_out, 'lat', attname = 'standard_name', attval = 'latitude', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'lat', attname = 'units', attval = 'degrees_north', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'lat', attname = '_CoordinateAxisType', attval = 'GeoY', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'lat', attname = 'axis', attval = 'Y', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'lon', attname = 'standard_name', attval = 'longitude', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'lon', attname = 'units', attval = 'degrees_east', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'lon', attname = '_CoordinateAxisType', attval = 'GeoX', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, 'lon', attname = 'axis', attval = 'X', prec = NA, verbose=FALSE, definemode=FALSE )

# Put global attributes in file
ncatt_put(nc_out, varid = 0, attname = 'title', attval = 'MCD43A Calculated Albedo based on RossThickLiSparseReciprocal BRDF model', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, varid = 0, attname = 'description', attval = 'Calculated from MCD43A1 BRDF/Albedo Model Parameters Product (MODIS/Terra BRDF/Albedo Model_1 Daily L3 Global 500m SIN Grid)\n
                                        https://www.umb.edu/spectralmass/terra_aqua_modis/v006/mcd43a1_brdif_albedo_model_parameters_product', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, varid = 0, attname = 'Conventions', attval = 'CF-1.6', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, varid = 0, attname = 'source1', attval = 'MCD43A1 NetCDF sourced from APPEEARS: https://lpdaac.usgs.gov/tools/data_access/appeears', prec = NA, verbose=FALSE, definemode=FALSE )
ncatt_put(nc_out, varid = 0, attname = 'source2', attval = 'Processed using albedo formula implemented in R; https://github.com/jjmcnelis/mcd43-calculated-albedo', prec = NA, verbose=FALSE, definemode=FALSE )

# Write out and close
nc_out <- nc_open(albedo_out, write = TRUE)
nc_close(nc_out)
