%% Response time (Tau) map
% This function will generate a response time map using data from the ouput
% workspace of "MonteCarlo_error_forLinearInversion.m". A shapefile may
% optionally be output for visualization in a GIS program like ArcMap or
% QGIS. This script requires the mapping toolbox of MATLAB.
%
% Required inputs:
%   1) Stau - response time of each streamnode in the river profile output
%             from "MonteCarlo_error_forLinearInversion.m"
%   2) S - STREAMobj of the stream network output from "MonteCarlo_error_forLinearInversion.m"
%   3) DEM - Digital elevation model or another GRIDobj.
%
% Optional inputs:
%   1) exportname - Name of exported shapefile. if exportname is included the 
%               function will output a shapefile of the stream network
%               with response time attributes to the working directory. 
%               The default option will not export a shapefile
%
% Outputs:
% MS - map object of the stream network with response time attributes.
% Shapefile - there is the option to export the map object to a shapefile.
%             The shapefile has polyline geometry and is aggregated to
%             segments of 3*DEM cellsize length with attributes for min and
%             max response time.
%
% Author: James A. Fisher [jamesfisher.gis@gmail.com]
% Last modified: 5/19/2021
%%
function [MS] = taumap(Stau,S,DEM,varargin)

p = inputParser;         
p.FunctionName = 'TauMap';

% required inputs
addRequired(p,'Stau', @(x) ismatrix(x));
addRequired(p,'S', @(x) isa(x,'STREAMobj'));
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));

addOptional(p, 'exportname', []);

parse(p,Stau, S, DEM, varargin{:});
Stau = p.Results.Stau;
S = p.Results.S;
DEM   = p.Results.DEM;
exportname = p.Results.exportname;

figure(1)
imageschs(DEM,[],'colorbar',false,'colormap',[.9 .9 .9]); hold on
scatter(S.x,S.y,2,Stau);
cb=colorbar;
cb.Label.String = 'Response Time (yrs)';

if isempty(p.Results.exportname)
    MS = STREAMobj2mapstruct(S,'seglength',S.cellsize*3,'attributes',{'mintau' Stau @min, 'maxtau' Stau @max});
else
    MS = STREAMobj2mapstruct(S,'seglength',S.cellsize*3,'attributes',{'mintau' Stau @min, 'maxtau' Stau @max});
    shapewrite(MS,exportname);
end
end
    

