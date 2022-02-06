%% ClipDEM2Drainage - Clip rasters to selected catchment boundary
% This script will clip GRIDobj rasters to a selected watershed boundary.
% Work through the script section by section.

% I recommend that you first clip the rasters to a box that encompasses
% the a single catchment rather than a regional GRIDobj to avoid long processing time.

% The input rasters must be in text format. Rasters can be exported to .txt
% format in ArcGIS using the tool "Raster to ASCII". If geology and
% area weighting grid are included, all rasters must be EXACTLY the same
% cellsize.

% Workflow:
% Load in .txt raster files as GRIDobjs.
% Calculate the river network.
% Using Segment picker select the outlet location of your catchment.
% Clip GRIDobjs to that catchment boundary
% If needed, manually edit the rasters to exactly the same size.
% Clipped GRIDobjs and a stream shapefile will be exported to the 
% working directory.

% If you are working with a DEM only, comment out the lines for geol and precip.

%Author: James A. Fisher; jaf319@lehigh.edu
%Last Modified: 4/27/2021

%% Load data
DEM = GRIDobj('name_demraw.txt');
geol = GRIDobj('name_geolraw.txt');
Aweight = GRIDobj('name_ppraw.txt'); %upstream area weighting grid. Could
%be annual precipitation or weighting grid to simulate area change.
%% Calculate stream network
FD  = FLOWobj(DEM,'preprocess','carve');
A   = flowacc(FD);
DA   = flowacc(FD);
S = STREAMobj(FD,'minarea',1e6/(FD.cellsize^2)); % 10 x 10 x 10000 = 1 km2
S = removeshortstreams(S,1e5);
%% Select catchment
% Press "Enter" then click on the drainage outlet. Only select 1 drainage at a time.
% When asked if you would like to continue picking streams click "No"
[Sc]=SegmentPicker(DEM,FD,A,S,1,'direction','up');
%% Clip GRIDobjs to catchment boundary
D = drainagebasins(FD,Sc);
figure(1)
imageschs(D); hold on
plot(S); title('catchment boundary'); hold off
if exist('DEM','var'); DEM = clip(DEM,D>0);
    figure(2)
    imageschs(DEM); hold on
    plot(Sc); title('clipped DEM'); hold off
end
if exist('geol','var'); geol = clip(geol,D>0);
    figure(3)
    imageschs(geol,[]); hold on
    plot(Sc); title('catchment geology'); hold off
end
if exist('Aweight','var'); Aweight = clip(Aweight,D>0);
    figure(4)
    imageschs(Aweight); hold on
    plot(Sc); title('Area weighting grid'); hold off 
end
%% Manually edit GRIDobjs the EXACTLY the same size
% This section may be needed if you are clipping geol and Aweight GRIDobjs.
% If you are only using a DEM skip this section.
% Double click to open each of the GRIDobjs that were just clipped. 
% The "Size" property of each GRIDobj must be EXACTLY the same. Sometimes
% there may be 1 or 2 row or column discrepancy of NaN values.

% The lines below are some commands to manually remove extra rows from clipped
% GRIDobjs to ensure that they are exactly the same size.
% Be careful, there may be a unique solution to each set of data.

%DEM.Z(:,end) = [];  %remove the last column from DEM
%DEM.Z(end,:) = [];  %remove the last row from DEM
%DEM.size(1,1) = length(DEM.Z(:,1)); %update size property - row
%DEM.size(1,2) = length(DEM.Z(1,:)); % column

%geol.Z(:,end) = []; %remove the last column from geol
%geol.Z(end,:) = []; %remove the last row from geol
%geol.size(1,1) = length(geol.Z(:,1)); % update size property - row
%geol.size(1,2) = length(geol.Z(1,:)); % column

%Aweight.Z(:,end) = []; %remove the last column from Aweight
%Aweight.Z(end,:) = []; %remove the last row from Aweight
%Aweight.size(1,1) = length(Aweight.Z(:,1)); %update size property -row
%Aweight.size(1,2) = length(Aweight.Z(1,:)); % column
%% Export clipped DEM, geology, precipiation rasters and stream shapefile
%Edit the name but keep the naming convention "name_ws.txt",
%"name_geol.txt", "name_Aweight.txt"
GRIDobj2ascii(DEM, 'name_ws.txt');
GRIDobj2ascii(geol, 'name_geol.txt');
GRIDobj2ascii(Aweight, 'name_Aweight.txt');
MS = STREAMobj2mapstruct(Sc,'seglength',100,'attributes',...
    {'uparea' (DA.*(DA.cellsize^2))});
shapewrite(MS,'name_streams.shp');