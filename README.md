# LongProfileLinearInversion-FisherEtal-2022

## Suite of MatLab scripts for linear inversion of fluvial topography with variable erodibility

This package contains scripts that can be used to process raster data and perform linear
inversion of fluvial topography following the methods put forth in Goren (2014) and Gallen (2018). This document outlines each script’s purpose and provides
detailed documentation of data structures used in the linear inversion Monte Carlo simulation.
An example dataset is included for the Enza catchment and the proper input parameters are
already set in `MainLinearInversion_MonteCarlo.m`


To run these scripts, you will need:
- An install of [Topotoolbox 2](https://topotoolbox.wordpress.com/) and MATLAB >= R2019b
- Rasters clipped to a box that encompasses the entire desired catchment in text file format. If
you are using ArcGIS, rasters can be exported to text files with the tool “Raster to ASCII”.
- An estimate for basin-averaged erosion rate.

### Description

#### `ClipDEM2Drainage.m`:
- This script will process rasters to the correct format to be input into the linear inversion
scheme. The rasters must be exactly the same size with the same cellsize. The processed
rasters will be clipped to the selected catchment boundary and be output as new text files to
be input into the linear inversion scheme.
- Processed raster data for each of the six catchments from this study are included.

#### `Get_Ksndistribution_wErrors.m`:
- This script is only necessary if comparing uplift histories from the linear inversion across
multiple catchments. The inversion must be run using consistent K values for each rock-type.
Running this script will calculate the least-squares ksn for each rock-type and output a table of
the distribution of ksn following the format [numerical ID, ksn, ksn undertainty, number of
datapoints]. These distributions are used to calculate the average regional K value for each
rock-type weighted by the fraction of data points.

#### `MainLinearInversion_MonteCarlo.m`:
This script is designed to run a linear inversion of fluvial topography for a specified
number of iterations using randomly sampled values of ksn (normalized steepness) and E (Basinwide erosion rate) within a uniform distribution of the error (+ or – calculated errors). If
comparing between multiple basins, one may input K values and errors for each rock type.

This script contains all the required functions and is configured to run with the example
data set from the Enza catchment. Run this program from a working directory that includes only
this script, the example data, and the [Topotoolbox 2](https://topotoolbox.wordpress.com/) MatLab library (Shwanghart & Scherler, 2014). See the
`linear_inversion_MonteCarlo` function header for more information.

##### Required Inputs:
- DEM -- digital elevation model clipped to the watershed boundary
- E –- basin-wide erosion rate in m/yr
- E_error –- basin-wide erosion rate uncertainty in m/yr
- n –- number of model iterations (set to 1 to run without error)

##### Optional Inputs:
- tau_inc – timestep increment for final uplift history. Default is 250kyr
- mn – concavity index (theta). Default is 0.5
- geol – geology raster with integer ID for each rock-type. The raster must be clipped to the
watershed boundary and must be exactly the same size and cellsize as the DEM.
 - Note: if there is no input geol raster the program will assume uniform rock-type.
- K – vertical 1 x q array of K values for each lithology (q = # of rock-types present)
 - Note: If you input K values the program will use the input values instead of calcuated K
values from ksn and E.
- K_error – vertical array of uncertainty associated with each input K value.
- Aweight – raster for weighting the upstream area parameter (a proxy for discharge). There is the
option to weight upstream area by a precipitation gradient or to simulate drainage are
loss/gain by weighting the headwater drainage area. The raster must be clipped to the
catchment boundary and must be exactly the same size and cellsize as the DEM and geol.
- meanA –This is used to calculate the drainage area weighting factor. In the case of an annual
precipitation raster, this would be the average precipitation rate. If weighting factors are
already calculated for the Aweight grid, set meanA = 1. If no meanA is given, the
program will use the average value in the raster.

##### Relevant Parameters:
- n – number of model runs (if n=1 the model will run once without including error)
- m – number of stream nodes.
- q – number of rock-type bins
- En – n x 1 array of randomly sampled erosion rates from withing a uniform distribution within +- the E_error.
- ksn_n – n x q matrix of randomly sampled ksn values for ‘q’ rock-types within the associated error windows.
- Kn – n x q matrix of random K values; Kn = En/Ksn_n
timesteps – 50 x n matrix that holds the output timesteps from ‘n’ inversion iterations in Myr. Cells not filled by time data are set to NaN. Uplift histories are limited to 50 timesteps.
- U_rates – 50 x n matrix that hold the output uplift rates in from ‘n’ iterations in mm/yr. Cells not filled by time data are set to NaN.
- Mean_Urate – The average U_rates at each timestep after ‘n’ iterations.
- Med_Urate – The median U_rates at each timestep after ’n’ iterations.
- BF_timesteps – timesteps to best-fit data. Array from 0 to the longest record from ‘n’ iterations.

DATA STRUCTURE IMAGES HERE

##### Workflow:
1. Initialize data and generate optional weighted flow accumulation.
2. Bin stream nodes by geology type
a. Geol is an optional input. If there is no geol called, the program will bin
all nodes into one bin – uniform rock type.
3. Run “invert_for_ksn” function to generate the forward model “A” matrix of Chi
values and least-squares ksn for each rock type (KsnMod) with errors (KsnStd).
a. There are optional inputs for K and K_error to use specified values for
erodibility rather than K calculated from least-squares ksn and E values.
4. Generate En and ksn_n matrices of random erosion and ksn values within the
respective error windows for ‘n’ iterations. If there are input K and K_errors,
generate a list of random K values (Kn). The first set in the list will be K values
calculated with the original ksn and E.
5. Run a for loop for ‘n’ iterations. Each iteration calculates Stau from Kn and the
forward ‘Achi’ chi matrix (Stau = Achi/Kn). Stau is input into.
“linear_inversion_block_uplift_variableK”. Output uplift and time data for each
iteration are saved in timesteps and U_rates.
6. Calculate Mean_Urate (mean) and Med_Urate (median) uplift rates at each
timestep.
7. Calculate the median length of uplift records from final_ts.
8. Plot best-fit data with uplift histories from ‘n’ iterations.

 ##### Example output graph from the Enza basin with n=100.
 
PLACE IMAGE HERE

#### `taumap.m`:
- This is an additional function that is designed to be run using the output workspace from `MainLinearInversion_MonteCarlo.m`. This function will generate a map object of the stream network in map view colored by response time. There is the option to export a shapefile of the stream network with response time attributes. This function can be called with examples from the `TauMap_U2T_Example.m` script.

#### `uplift2topography.m`:
- This is an additional function that is designed to be run using the output workspace from `MainLinearInversion_MonteCarlo.m`. The function will integrate uplift with respect to response time to model long profile elevation for each timestep. The function will also plot the stream network in each timestep in map view. There is an optional setting to output .gif animations of the long profile and map view evolution. This function can be called with examples from the `TauMap_U2T_Example.m` script.

### References Cited:
Gallen, S. F. (2018). Lithologic controls on landscape dynamics and aquatic species
<br>&nbsp;&nbsp;&nbsp;&nbsp; evolution in post-orogenic mountains. *Earth and Planetary Science Letters, 
<br>&nbsp;&nbsp;&nbsp;&nbsp; 493, 150–160. https://doi.org/10.1016/j.epsl.2018.04.029*

Goren, L., Fox, M., & Willett, S. D. (2014). Tectonics from fluvial topography using
<br>&nbsp;&nbsp;&nbsp;&nbsp; formal linear inversion: Theory and applications to the Inyo Mountains, California.
<br>&nbsp;&nbsp;&nbsp;&nbsp; Journal of *Geophysical Research: Earth Surface, 119(8), 1651–1681. https://doi.org/10.1002/2014JF003079*

Schwanghart, W., & Scherler, D. (2014). Schwanghart, W., & Scherler, D. (2014).
<br>&nbsp;&nbsp;&nbsp;&nbsp; Short Communication: TopoToolbox 2 – MATLAB-based software for topographic 
<br>&nbsp;&nbsp;&nbsp;&nbsp; analysis and modeling in Earth surface sciences. *Earth Surface Dynamics,
<br>&nbsp;&nbsp;&nbsp;&nbsp; 2(1), 1–7. Copernicus GmbH. Retrieved from https://doi.org/10.5194%2Fesurf-2-1-2014*
