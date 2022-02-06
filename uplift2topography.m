%% Integrate uplift history for topography 
% this script will integrate uplift history with respect to response time
% for each timstep following the methods of Goren et al (2014). This
% function is designed to be run from the output workspace of
% "MonteCarlo_error_forLinearInversion.m".
%
% Required inputs:
%   1) Uhist - uplift history, either Med_Urate or Mean_Urate
%   2) timesteps - the timesteps associated with the uplift history,
%                 BF_timsteps
%   3) Stau - response time of each streamnode in the river profile output
%             from "MonteCarlo_error_forLinearInversion.m"
%   4) S - STREAMobj of the stream network output from "MonteCarlo_error_forLinearInversion.m"
%   5) DEM - elevation GRIDobj
%
% Optional inputs:
%   1) export - input 'true' to export .gif animations.
%
% Output:
%   1) z_t - modeled present day elevation.
%   2) optional .gifs - ouput to the working direct with the name
%   profile.gif and map.gif.
%
% Author: James A. Fisher - (jamesfisher.gis@gmail.com)
% Last modified: 5/24/2021
function [z_t] = uplift2topography(Uhist,timesteps,Stau,S,DEM,varargin)

p = inputParser;         
p.FunctionName = 'uplift2topography';

% required inputs
addRequired(p,'Uhist', @(x) ismatrix(x));
addRequired(p,'timesteps', @(x) ismatrix(x));
addRequired(p,'Stau', @(x) ismatrix(x));
addRequired(p,'S', @(x) isa(x,'STREAMobj'));
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));

addOptional(p, 'export', []);

parse(p,Uhist, timesteps, Stau, S, DEM, varargin);
Uhist = p.Results.Uhist;
timesteps = p.Results.timesteps;
Stau = p.Results.Stau;
S = p.Results.S;
DEM   = p.Results.DEM;
export = p.Results.export;

%Prep data for integration
label = flip(string(timesteps));
TSflip = flip(timesteps);
timesteps(1,:) = 0.0001; %set first timestep to 100yrs to make integration work
Uhist = flip(Uhist/1e3);
Sz = DEM.Z(S.IXgrid);
stream = [S.x,S.y,S.distance,(Sz-min(Sz)),Stau];
stream = sortrows(stream,5); % sort according the increasing Stau
cols = prism(length(timesteps));

%Prep figures
h = figure(1); hold on
m = figure(2); hold on
figure(1)
subplot(2,1,1);
plot(stream(:,3),stream(:,4),'.','Color',[0.8 0.8 0.8]); hold on 
figure(2)
imageschs(DEM,[],'colorbar',false,'colormap',[.9 .9 .9]); hold on

%Select subset of data for each integration step and interate for each
%timestep
for i = 2:length(timesteps)
    U_rate = Uhist(i);
    tau_lim = timesteps(i)*1e6;
    ind = stream(:,5) <= tau_lim; %timestep index
    stream_t = stream(ind,:);
    Stau_t = stream_t(:,5);
    U_t = [];
    for t = 1:length(Stau_t) % for loop creating an array with constant U for each node
        %reminder that this is assuming block uplift. For a flexural response,
        %a diffusion function with distance could be used for each time step
        U_t = [U_t;U_rate];
    end
    
    % Integrate uplift with respect to response time
    z_t = cumtrapz(Stau_t,U_t); 
    z_t = z_t + min(stream(:,4));
    
    %Plot modeled long profile for each timestep 
    figure(1)
    subplot(2,1,1);
    plot(stream_t(:,3),z_t,'.','MarkerSize',1,'Color',cols(i,:));
    legend({'Present day'});
    xlabel('Distance (m)'); ylabel('Elevation (m)');
    title(strcat('Modeled elevation: ',label{i},' Myr BP'));
    subplot(2,1,2);
    stairs(TSflip*1e6,Uhist,':*b','lineWidth',2); hold on
    xlabel('Myr BP'); ylabel('U (m/yr)');
    plot(TSflip(i)*1e6,Uhist(i),'-o','MarkerSize',15,'Color',cols(i,:));
    
    %Plot the data in map view
    figure(2)
    L = plot(stream_t(:,1),stream_t(:,2),'.','Color',cols(i,:)); 
    title(strcat('knickpoint propagation: ',label{i},' Myr BP'));
    if i > 1
        uistack(L,'down',i-2); 
    end
    
% if record is chosen, make gif animation in map and long profile view and
% export
    if isempty(p.Results.export)
        pause(1);
    else
        pause(1);
        frame = getframe(h); %profile
        framem = getframe(m); %map
        im = frame2im(frame); %prfile
        imm = frame2im(framem); %map
        [imind,cm] = rgb2ind(im,256); %profile
        [imindm,cmm] = rgb2ind(imm,256); %map
  %export .gif animation of the profile
        if i == 2
            imwrite(imind,cm,'profile.gif','Loopcount',inf,'DelayTime',2); %profile
            imwrite(imindm,cmm,'map.gif','Loopcount',inf,'DelayTime',2); %map
        else 
            imwrite(imind,cm,'profile.gif','WriteMode','append','DelayTime',2); %profile
            imwrite(imindm,cmm,'map.gif','WriteMode','append','DelayTime',2); %map
        end 
    end
end
end