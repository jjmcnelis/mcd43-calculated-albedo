#####################################################################################
# Function for looping through albedo parameters and timesteps
#####################################################################################

proc_albedo <- function(band){
  
  # Initialize output albedo stacks
  black <- stack()
  white <- stack()
  actual <- stack()
  
  if(band == 'vis'){
    lyr <- 'BRDF_Albedo_Parameters_vis'
    lyrqc <- 'BRDF_Albedo_Band_Mandatory_Quality_vis'

  }
  
  if(band == 'nir'){
    # Get NIR variable
    lyr <- 'BRDF_Albedo_Parameters_nir'
    # Get NIR QC variable
    lyrqc <- 'BRDF_Albedo_Band_Mandatory_Quality_nir'
  }
  
  if(band == 'shortwave'){
    # Get shortwave variable
    lyr <- 'BRDF_Albedo_Parameters_shortwave'
    # Get shortwave QC variable
    lyrqc <- 'BRDF_Albedo_Band_Mandatory_Quality_shortwave'
  }
  
  # Get variable
  var <- nc$var[[lyr]]
  # Get QC variable
  varqc <- nc$var[[lyrqc]]
  
  # Get _FillValue
  fillvalue <- ncatt_get(nc, lyr, "_FillValue")
  
  # Get variable size, dimension count (4), and timestep count
  varsize <- var$varsize
  varndims <- var$ndims
  
  # Get size and dimension count
  qcndims <- varqc$ndims
  qcsize <- varqc$varsize
  
  # Loop through timesteps
  for(tstep in 1:nt){
    print(paste("Processing", band, "timestep", tstep, "of", nt, ". . .", sep = " "))
    
    # Convert tstep to real date
    dt <- as.Date(nc_time[tstep], origin = as.Date("2000-01-01"))
    # Calculate array of solar zenith angles for each pixel for timestep n
    # I NEED TO FIGURE OUT WHAT THE DAAC AND MCD43 DO FOR THIS!! WHY ARE ALL VALUES BELOW 90? WHY DOES IT ONLY TAKE LATITUDE? THE LUTEXE??
    SZN <- (mapply(SZA, timein = rep(dt, nrow(nc_latlon)), Lat = nc_latlon[,1], Lon = nc_latlon[,2]) - 90)
    # Get array of LUT values for each SZN
    LUTarr <- as.vector(mapply(LUT, band = rep('vis', length(SZN)), szn = SZN, sod = rep(SolOptDepth, length(SZN))))
    
    # Initialize start and count to read one timestep of the variable
    vstart <- rep(1,varndims)
    vstart[varndims] <- tstep
    vcount <- varsize
    vcount[varndims] <- 1
    
    qstart <- rep(1,qcndims)
    qstart[qcndims] <- tstep
    qcount <- qcsize
    qcount[qcndims] <- 1
    
    # Read visible data for timestep n from nc structure
    vardata <- ncvar_get(nc, var, start=vstart, count=vcount)
    # Read visible qc for timestep n from nc structure
    qcdata <- ncvar_get(nc, varqc, start=qstart, count=qcount)
    
    # Set _FillValue to NA
    vardata[vardata == fillvalue$value] <- NA
    
    # Get 3 parameters as raster objects
    vardata1 <- raster(t(vardata[1,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
    vardata2 <- raster(t(vardata[2,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
    vardata3 <- raster(t(vardata[3,,]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
    
    # Apply QC mask to vardata (valid == 0 or 1)
    vardata[!qcdata %in% c(0,1)] <- NA
    
    # Calculate Black Sky Albedo for timestep n
    nBlack <- BlackSA(vardata1, vardata2, vardata3, SZN)
    
    # Calculate White Sky Albedo for timestep n
    nWhite <- WhiteSA(vardata1, vardata2, vardata3)

    # Calculate Actual Albedo for timestep n
    nActual <- ActualSA(nWhite, nBlack, LUTarr)
    
    # Add albedos to raster stacks
    black <- stack(black, nBlack)
    white <- stack(white, nWhite)
    actual <- stack(actual, nActual)
    
  }
  
  return(c(black, white, actual))
  
  # Remove variables
  remove(var, varqc, fillvalue, varsize, varndims, qcndims, qcsize, dt, SZN, LUTarr, vstart, vcount, 
         qstart, qcount, vardata, qcdata, vardata1, vardata2, vardata3, nBlack, nWhite, nActual)
  
}
