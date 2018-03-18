# R Scripts for Calculating Black, White, and Actual Albedo (MCD43A) from MCD43A1

### Why?

While some users are content to use the MODIS black-sky albedo at local solar noon and the white-sky albedo measures as provided in [**MCD43A3**](https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd43a3_v006), most researchers want to make use of the BRDF model parameters ([**MCD43A1**](https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd43a1_v006)) to obtain black-sky albedos at other solar illumination angles through the use of a simple polynomial or to combine the black-sky and white-sky albedos as a function of optical depth to calculate "actual" or blue-sky albedos for validation studies.

[**MCD43A**](https://modis.ornl.gov/documentation.html#MCD43) is the unofficial MODIS product code for **black-sky**, **white-sky**, and **actual (blue-sky) albedo** calculated using the BRDF model parameters and user-specified **solar optical depth** and **solar zenith angle**. The MCD43A Calculated Albedo Product is available through the [MODIS Global Subsetting and Visualization Tool](https://modis.ornl.gov/cgi-bin/MODIS/global/subset.pl) maintained by the [ORNL DAAC](https://daac.ornl.gov). Alternatively, it can be calculated from MODIS HDF files using the binaries available at the [product page](https://www.umb.edu/spectralmass/terra_aqua_modis) maintained by the PI of MCD43, Dr. Crystal Schaaf.

In MODIS Collection 5, MCD43 was generated in 16-day composite time intervals. In Collection 6 (C6), MCD43 is a daily product; thus, the data volume is prohibitively large. The large data volume for MCD43 in C6 was the motivation for writing these scripts. 

See below for more information about MCD43A Calculated Albedo.

For more information about MODIS, visit the LP DAAC's [MODIS Overview](https://lpdaac.usgs.gov/dataset_discovery/modis) page.

### How?

#### Input data from AppEEARS

The albedo formulas are implemented in R. They take as input a subset of the MCD43A1 product in NetCDF format generated through [AppEEARS](https://lpdaac.usgs.gov/tools/data_access/appeears). AppEEARS is a data delivery tool maintained by the LP DAAC offering spatial and temporal subsets of MODIS data and other data products. You will need to sign up for a NASA Earthdata account to use AppEEARS.

To request a compatible subset, choose **Area Sample** from the ***Extract*** button in AppEEARS. Click **Start a new request**. 
1. Provide a name for your sample.
2. Choose a spatial extent by either uploading an ESRI Shapefile or GeoJSON, or by drawing a box/polygon on the map.
3. Choose a temporal extent for the subset.
4. Search for MCD43A1 under ***Select the layers to include in the sample*** and choose the C6 product: **MCD43A1.006**
5. Select the following layers:
	*	**BRDF_Albedo_Parameters_vis**
	*	**BRDF_Albedo_Parameters_nir**
	*	**BRDF_Albedo_Parameters_shortwave**
	*	**BRDF_Albedo_Band Mandatory_Quality_vis**
	*	**BRDF_Albedo_Band Mandatory_Quality_nir**
	*	**BRDF_Albedo_Band Mandatory_Quality_shortwave**
6. For ***File Format***, choose **NetCDF**.
7. For ***Projection***, choose **Geographic**.

You will receive an email upon order completion.

An example input file is provided: *./example-data/MCD43A1.006_aid0001.nc*. You can submit an identical AppEEARS request by uploading the supplied GeoJSON: *./example-data/MCD43A1-FL-Everglades-request.json*

#### Usage

The primary script (*main_albedo.R*) calls the functions defined in *f_albedo.R* and *proc_albedo.R*. Open *main_albedo.R* and make a few edits to the settings in the first section:
```
# Set your working directory
# Contains R scripts, lookup tables (vis_lut.csv, nir_lut, sw_lut.csv), and input NetCDF
wkd <- '<<main_albedo.R file location>>'
setwd(wkd)

# Set path and filename for input NetCDF. Default output from AppEEARS should be MCD43A1.006_aid0001.nc
nc_in = paste(wkd, 'example-data/', 'MCD43A1.006_aid0001.nc', sep = '')

# Set solar optical depth (value between 0.00 and 0.98 x0.02 with 2-digit precision)
SolOptDepth <- 0.20
# Set solar zenith angle (0.0-89.0 degrees, or "local" for local solar noon)
SolZenAngle <- 'local'

# Set path to output NetCDF
out_dir <- paste(wkd, 'output/', sep = '')
out_filename <- ''
suppressWarnings(dir.create(out_dir))
```

Execute the script via command line:
```
> Rscript main_albedo.R

In Windows, you may need to point to your Rscript executable:
> C:\Path\to\R\directory\Rscript.exe main_albedo.R
```
Output will be a NetCDF with the following data variables:
* **visible_black_sky_albedo**
* **visible_white_sky_albedo**
* **visible_actual_albedo**
* **nir_black_sky_albedo**
* **nir_white_sky_albedo**
* **nir_actual_albedo**
* **shortwave_black_sky_albedo**
* **shortwave_white_sky_albedo**
* **shortwave_actual_albedo**

An example output file is provided: *./example-data/output/MCD43A_calcalbedo.nc*

## MCD43A Calculated Albedo

The primary resource for information about MCD43 products is the [product page](https://www.umb.edu/spectralmass/terra_aqua_modis/v006) maintained by Dr. Schaaf's group at UMass Boston. As the only institution offering on-demand delivery of MCD43A Calculated Albedo, the ORNL DAAC has a [web page](https://modis.ornl.gov/documentation.html#MCD43) with a more concise explanation of its derivation.

Albedo is defined as the ratio of upwelling to downwelling radiative flux at the surface. Downwelling flux may be written as the sum of a direct component and a diffuse component. Black-sky albedo (directional hemispherical reflectance) is defined as albedo in the absence of a diffuse component and is a function of solar zenith angle. White-sky albedo (bihemispherical reflectance) is defined as albedo in the absence of a direct component when the diffuse component is isotropic. Black-sky albedo and white-sky albedo mark the extreme cases of completely direct and completely diffuse illumination. Actual albedo is a value which is interpolated between these two as a function of the fraction of diffuse skylight which is itself a function of the aerosol optical depth (Lewis and Barnsley 1994; Román et al. 2010; Schaaf et al. 2002).

The MCD43A1 BRDF/Albedo Model Parameters Product (MODIS/Terra BRDF/Albedo Model_1 Daily L3 Global 500m SIN Grid) supplies the weighting parameters associated with the RossThickLiSparseReciprocal BRDF model that best describes the anisotropy of each pixel. These three parameters ( fiso , fvol , fgeo ) are provided for each of the MODIS spectral bands as well as for three broad bands (0.3-0.7µm, 0.7-5.0µm, and 0.3-5.0µm). These parameters can be used in a forward version of the model to reconstruct the surface anisotropic effects and thus correct directional reflectances to a common view geometry (this is the procedure that is used to produce MCD43A4 Nadir BRDF-Adjusted Reflectances - NBAR) or to compute the integrated black-sky (at some solar zenith angle) and white-sky albedos (as are done for MCD43A3). Alternately, the parameters can be used with a simple polynomial to easily estimate the black-sky albedo with good accuracy for any desired solar zenith angle. The polynomial is as follows:

![Broken Image Link --- black sky albedo](https://github.com/jjmcnelis/mcd43-calculated-albedo/blob/master/readme/blackskyalbedo.jpg?raw=true)

The appropriate constants are:

| Term   | Isotropic (iso) | RossThick (vol) | LiSparseR (geo) |
|--------|-----------------|-----------------|-----------------|
| *g0*   | 1.0             | -0.007574       | -1.284909       |
| *g1*   | 0.0             | -0.070987       | -0.166314       |
| *g2*   | 0.0             | 0.307588        |  0.041840       |

Similarly, the white-sky albedo can be computed by using the equation:

![Broken Image Link --- white sky albedo](https://github.com/jjmcnelis/mcd43-calculated-albedo/blob/master/readme/whiteskyalbedo.jpg?raw=true)

and the estimates of the white-sky kernel integrals are:

| Term                     | Isotropic (iso) | RossThick (vol) | LiSparseR (geo) |
|--------------------------|-----------------|-----------------|-----------------|
| White-sky integral *g*   | 1.0             | 0.189184        | -1.377622       |

Actual (blue-sky) albedo is interpolated between black- and white-sky albedo as a function of the fraction of diffuse skylight which is itself a function of the aerosol optical depth:

![Broken Image Link --- actual albedo](https://github.com/jjmcnelis/mcd43-calculated-albedo/blob/master/readme/actualalbedo.png?raw=true)

Lewis, P., and M. J. Barnsley, Influence of the sky radiance distribution on various formulations of the earth surface albedo, in Proc. Conf. Phys. Meas. Sign. Remote Sens., Val d'Isere, France, pp. 707-715, 1994.

Román, M. O., C. B. Schaaf, P. Lewis, F. Gao, G. P. Anderson, J. L. Privette, A. H. Strahler, C. E. Woodcock, M. Barnsley, Assessing the coupling between surface albedo derived from MODIS and the fraction of diffuse skylight over spatially-characterized landscapes, Remote Sensing of Environment, 114, 738-760,2010.

Schaaf, C. B., F. Gao, A. H. Strahler, W. Lucht, X. Li, T. Tsang, N. C. Strugnell, X. Zhang, Y. Jin, J.-P. Muller, P. Lewis, M. Barnsley, P. Hobson, M. Disney, G. Roberts, M. Dunderdale, C. Doll, R. d'Entremont, B. Hu, S. Liang, and J. L. Privette, First Operational BRDF, Albedo and Nadir Reflectance Products from MODIS, Remote Sens. Environ., 83, 135-148, 2002.
