#####################################################################################
# These functions are called by the loop function in main_albedo.R
#
#::METHODS::
# 
#   https://modis.ornl.gov/documentation.html
#   https://www.umb.edu/spectralmass/terra_aqua_modis/modis_brdf_albedo_product_mcd43
#
# Black-sky Albedo = 
#   Parameters_01 + 
#   Parameters_02 * (-0.007574 + (-0.070987 * szn^2) + (0.307588 * szn^3)) + 
#   Parameters_03 * (-1.284909 + (-0.166314 * szn^2) + (0.041840 * szn^3))
#
# White-sky Albedo = 
#   Parameters_01 + 
#   Parameters_02 * (0.189184) + 
#   Parameters_03 * (-1.377622) 
#
# Actual (Blue-sky) Albedo = 
#   White-sky Albedo * f(optical depth, solar zenith angle, aerosol type, band) + 
#   Black-sky Albedo * (1 - f(optical depth, solar zenith angle, aerosol type, band)) 
#
#####################################################################################
# Albedo calculation constants
#####################################################################################

## Black sky constants

# Isotropic constant
g0iso <- 1.0
g1iso <- 0.0
g2iso <- 0.0

# RossThick constant
g0vol <- -0.007574
g1vol <- -0.070987
g2vol <- 0.307588

# LiSparseR constant
g0geo <- -1.284909
g1geo <- -0.166314
g2geo <- 0.041840

## White sky constants

gIso <- 1.0       # Isotropic
gVol <- 0.189184  # RossThick
gGeo <- -1.377622 # LiSparseR

#####################################################################################
# Albedo functions
#####################################################################################

# Convert SZN degrees to radians
Deg2Rad = 3.1415926535/180
SF = 1 # Scale factor (0.001) already applied to values by ncdf4. Set to 1.

# Black-sky albedo formula
BlackSA <- function(p1arr, p2arr, p3arr, szn){
  sznrad = szn*Deg2Rad
  return((p1arr*SF)+(p2arr*SF)*(g0vol+(g1vol*sznrad^2)+(g2vol*sznrad^3))+(p3arr*SF)*(g0geo+(g1geo*sznrad^2)+(g2geo*sznrad^3)))
}

# White-sky albedo formula
WhiteSA <- function(p1arr, p2arr, p3arr){
  return((p1arr*SF)*gIso+(p2arr*SF)*gVol+(p3arr*SF)*+gGeo)
}

# Actual albedo formula
ActualSA <- function(WSA, BSA, LUTVal){
  return((WSA*LUTVal)+(BSA*(1-LUTVal)))
}

#####################################################################################
# Lookup tables
#####################################################################################

LUT <- function(band, szn, sod){
  if(band == 'vis'){
    return(vislut[paste('X', sprintf('%.2f', sod), sep = ''),round(abs(szn), digits = 0)+1])
  }
  if(band == 'nir'){
    return(nirlut[paste('X', sprintf('%.2f', sod), sep = ''),round(abs(szn), digits = 0)+1])
  }
  if(band == 'shortwave'){
    return(swlut[paste('X', sprintf('%.2f', sod), sep = ''),round(abs(szn), digits = 0)+1])
  }
}
