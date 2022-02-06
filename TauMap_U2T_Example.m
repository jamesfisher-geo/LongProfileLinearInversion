%% Example Script for calling functions to generate a response time map and integrate for topography
% The two functions called in this script are designed to be run with data
% output from the "MonteCarlo_error_LinearInversion.m" script. 

% Author: James A. Fisher - jamesfisher.gis@gmail.com
%% Response time map - option to export shapefile
FD = FLOWobj(DEM);
S  = STREAMobj(FD,'minarea',1e6/(DEM.cellsize)^2);
exportname = 'StironeResponseTime'; %name for the output shapefile
%[MS] = taumap(Stau,S,DEM);
[MS] = taumap(Stau,S,DEM,'exportname', exportname); %export .shp to
%directory
%% generate modeled long profile and fluvial evolution by integrating uplift rate w.r.t response time
%[z_t] = uplift2topography(Med_Urate, BF_timesteps, Stau, S, DEM);
[z_t] = uplift2topography(Med_Urate, BF_timesteps, Stau, S, DEM,'export','true'); %will export .gif animations of long profile and map view