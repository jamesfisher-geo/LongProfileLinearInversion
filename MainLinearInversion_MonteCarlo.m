%% Monte Carlo style error propagation for linear inversion
% This script will run the linear_inversion_ksn_and_tau_variableK.m
% function that will run a linear inversion with variable erodibility by
% rock type.

% See the description documentation for outline of the workflow and data
% structure.
%
%Author:James A. Fisher jamesfisher.gis@gmail.com
%Sub-functions by: Sean F. Gallen; sean.gallen@colostate.edu
%
%Last Modified: 4/27/21
%% Initialize parameters
n = 1; %number of iterations. If you do not want to run the monte carlo, set n=1
% n=1 will run the model using the input values of E and K
E = 0.000519; % (m/yr) basin wide erosion rate 
E_error = 0.000129; % + or - error of erosion rate;
mn = 0.45; %Theta concavity value (m/n) 
K =  [1.88E-05;	6.67E-06;	7.31E-06;	1.08E-05;	1.48E-05;	2.90E-05;	4.86E-06]; %optional input
K_error = [1.85E-06;	8.05E-07;	7.64E-07;	1.31E-06;	1.09E-06;	8.28E-06;	1.06E-06]; %optional input
%meanA = 1; %Used to calculate upstream area weighting factors. Optional 
% name_ws.txt (DEM), name_geol.txt (geol), name_Aweight.txt (Area weighting grid).
wsname = 'enza'; % ENTER BASIN NAME - the input and outputs will auto-fill 
DEM = GRIDobj(strcat(wsname,'_ws.txt'));
geol = GRIDobj(strcat(wsname,'_geol.txt')); %Optional
%Aweight = GRIDobj(strcat(wsname,'_Aweight.txt')); weighting grid for upstream area. Optional
DEM = fillsinks(DEM);
%% Run the Monte Carlo - Run this section
% with uniform geology
%[timesteps,U_rates,Kn,Mean_Urate,Med_Urate,BF_timesteps,KsnMod,Stau,Atau,Schi,Achi,S] = linear_inversion_MonteCarlo(DEM,E,E_error,n,'mn',mn);
% with input geol
[timesteps,U_rates,Kn,Mean_Urate,Med_Urate,BF_timesteps,KsnMod,Stau,Atau,Schi,Achi,S] = linear_inversion_MonteCarlo(DEM,E,E_error,n,'mn',mn,'geol',geol);
% with input geol and K values
%[timesteps,U_rates,Kn,Mean_Urate,Med_Urate,BF_timesteps,KsnMod,Stau,Atau,Schi,Achi,S] = linear_inversion_MonteCarlo(DEM,E,E_error,n,'mn',mn,'geol',geol,'K',K,'K_error',K_error);
% With input geol, K, and Aweight 
%[timesteps,U_rates,Kn,Mean_Urate,Med_Urate,BF_timesteps,KsnMod,Stau,Atau,Schi,Achi,S] = linear_inversion_MonteCarlo(DEM,E,E_error,n,'mn',mn,'geol',geol,'K',K,'K_error',K_error,'Aweight',Aweight,'meanA',meanA);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Monte carlo inversion function
function [timesteps,U_rates,Kn,Mean_Urate,Med_Urate,BF_timesteps,KsnMod,Stau,Atau,Schi,Achi,S] = linear_inversion_MonteCarlo(DEM,E,E_error,n,varargin)
% This script uses "invert_for_ksn.m" and
% linear_inversion_block_uplift_variableK.m" scripts to invert for ksn and 
% invert response time for uplift rate. 
% read the help menu under "invert_for_ksn" and 
% linear_inversion_block_uplift_variableK.m" before starting. 
%
% DEM, geol, and Aweight GRIDobj's must be clipped to the the watershed
% extent and be EXACTLY the same size (row x column). Use the
% "ClipDEM2Drainage.m" script to prepare GRIDobj and export .txt files
% ready for input.
%
%
% Required Inputs:
%
%       (1) DEM - clipped DEM to the watershed boundary
%       (2) E - basin wide erosion rate
%       (3) E_error - basin wide erosion rate uncertainty
%       (4) n - number of iterations
%       
% Optional Inputs:
%
%       (1) tau_inc - time step for tau inversion. 250,000yrs is recommended.
%       (2) mn - best-fit theta value for the whole stream network. See "best_mn_with_chi.m".
%       (3) geol - clipped geology raster to the watershed boundary with
%           a numerical identifier for each lithology.
%               Note: if there is no input geol, program will assume
%               uniform geology
%       (4) K - Vertical array of K values for each lithology. 
%               Note: if there is no input K, the program will calculate K
%               values.
%       (5) K_error - Vertical array of K errors for each lithology
%       (6) Aweight - GRIDobj for weighting the upstream area parameter.
%       (7) meanA - weighted mean value for calulating a weighting factors.
%           If you enter an Aweight GRIDobj but don't specify an meanA, 
%           the script will find the average of Aweight.Z.
%               Note: the Aweight GRIDObj acts as a weighting factor to
%               upstream area. This can be defined by a precipitation
%               gradient or could be used to simulate basin area
%               loss/gain. If the Aweight GRIDobj already has the desired
%               weighting factors assigned, keep meanA = 1.
%       (8) crita - critical upstream area threshold for channelized
%           rivers. Default is 1e6 (1km^2). 
% Outputs:
%
%       (1) timesteps - 50 x n matrix that holds the output timesteps from
%                   ‘n’ inversion iterations in Myr. Cells not filled by time data
%                    are set to NaN.
%       (2) U_rates - 50 x n matrix that hold the output uplift rates in 
%                   from ‘n’ iterations in mm/yr. Cells not filled by time data 
%                   are set to NaN
%       (3) Kn - n x q matrix of random K values used in each iteration.
%                   The first set in the list are the input K or the best-fit
%                   calculated K values.
%       (4) Mean_Urate - The average U_rates at each timestep after ‘n’ iterations
%       (5) Med_Urate - The median U_rates at each timestep after ’n’ iterations
%       (6) BF_timesteps - timesteps of best-fit data. Array from 0 to the
%                   median final timestep from ‘n’ iterations
%       (7) KsnMod - least-squares ksn. If K is input this is backcalculated
%                  from K = E/ksn.
%       (8) Stau - array of response time values for each stream node.
%       (9) Atau - forward model matrix of change in response time at each
%                   stream node
%       (10)Schi - array of Chi values for each stream node.
%       (11)Achi - forward model matrix of change in Chi at each stream
%                  node
%       (12)S - STREAMobj used for inversion.

% Workflow:
%  Calculate flow accumulation grid (option to use a Aweight weighting)
%  and upstream area. > Generate the S_discrete variable (a stream network 
%  array with values for underlying geology) > invert for ksn binned by 
%  geology > generate n random numbers within error windows for E and Ksn > 
%  Calculate array of nxq erodibility values (K) > iterate n times - calculate 
%  a network of tau values with variable K (Stau) > and invert Stau for 
%  uplift history (U_rates and timesteps).
%
%   **See readme file for detailed description of data structures**

%Author:James A. Fisher jamesfisher.gis@gmail.com
%Sub-functions by: Sean F. Gallen; sean.gallen@colostate.edu
%
%Last Modified: 4/27/21

%%Parse Inputs
p = inputParser;         
p.FunctionName = 'linear_inversion_MonteCarlo';

% required inputs

addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'E', @(x) isscalar(x));
addRequired(p,'E_error', @(x) isscalar(x));
addRequired(p,'n', @(x) isscalar(x)); %number of iterations


% optional inputs
addOptional(p,'tau_inc', 2.5e5, @(x) isscalar(x));
addOptional(p,'mn', 0.5, @(x) isscalar(x));
addOptional(p,'geol',[],@(x) isa(x,'GRIDobj'));
addOptional(p,'K',[], @(x) ismatrix(x));
addOptional(p,'K_error',[], @(x) ismatrix(x)); %array of errors for Kmod
addOptional(p,'Aweight',[], @(x) isa(x, 'GRIDobj'));
addOptional(p,'meanA',[],@(x) isscalar(x));
addOptional(p,'crita',1e6,@(x) isscalar(x));
%addOptional(p,'flowOption', []);

parse(p,DEM, E, E_error, n, varargin{:});

DEM   = p.Results.DEM;
E   = p.Results.E;
E_error   = p.Results.E_error;
n = p.Results.n;

tau_inc = p.Results.tau_inc;
mn = p.Results.mn;
geol   = p.Results.geol;
K = p.Results.K;
K_error = p.Results.K_error;
Aweight = p.Results.Aweight;
meanA = p.Results.meanA;
crita = p.Results.crita;
%% Prepare flow accumulation and upstream area data
if isempty(p.Results.Aweight)
    FD = FLOWobj(DEM);
    DA = flowacc(FD)*(FD.cellsize^2);
else
    if isempty(p.Results.meanA)
       meanA = nanmean(Aweight.Z,[1,2,3]);
       Aweight.Z = (Aweight.Z)*(meanA^-1); % Weighting factor as fraction of the mean
       FD = FLOWobj(DEM);
       Aweight.refmat = FD.refmat;
       DA = flowacc(FD,Aweight)*(FD.cellsize^2);
    else
       Aweight.Z = (Aweight.Z)*(meanA^-1); % Weighting factor as fraction of the mean
       FD = FLOWobj(DEM);
       Aweight.refmat = FD.refmat;
       DA = flowacc(FD,Aweight)*(FD.cellsize^2);
    end
end
S  = STREAMobj(FD,'minarea',crita/(FD.cellsize^2));
%S = trunk(S); %uncomment this to only invert only the trunk channel. Must
    %also uncomment the S = trunk(s)  on line 806 in "invert_for_ksn.m"
S = klargestconncomps(S,1);
Sz = DEM.Z(S.IXgrid);
Sz = Sz - min(Sz);
S_DA = DA.Z(S.IXgrid);
%% Binning by geology and calculate ksn by rock-type
if isempty(p.Results.geol)
    S_discrete = ones(size(S.IXgrid)); %Set uniform geology
else
    S_discrete = geol.Z(S.IXgrid); %or bin streamnodes by underlying rock ID
end

if isempty(p.Results.K)  
    g_bins = unique(S_discrete); %Bin values set to unique values in S_discrete
% run non-negative least-squares inversion and plot results 
    [KsnMod,KsnStd,Achi,Schi] = invert_for_ksn(S,Sz,S_DA,S_discrete,'mn',mn);
% plot the data
    figure(2)
    subplot(1,2,1)
    errorbar(g_bins,KsnMod,KsnStd,'ko'); hold on
    xlabel('Unit ID'); ylabel('k_{sn}');
    SSz = Achi*KsnMod;
    subplot(1,2,2)
    plot(S.distance,Sz,'.','color',[0.7 0.7 0.7]); hold on
    plot(S.distance,SSz,'k.');
else % run the same inversion but plot the modeled ksn and modeled topography
    g_bins = unique(S_discrete); 
% run the non-negative least-squares inversion and plot results
    [Achi,Schi] = MakeDiscreteChiAMatrix(S,S_DA,S_discrete,mn);
    if numel(K)>numel(g_bins)
        error('There are more K values than geology types');
    end
    if numel(K)<numel(g_bins)
        error('There are less K values than geology types');
    end
    KsnMod = E./K; % temporary variable of the Ksn values from input K
    % propagate errors through ksn = E/K
    ksnerrorinputs = sqrt((E_error./E).^2+(K_error./K).^2).*KsnMod; % temporary variable of the Ksn errors from input K
    figure(2)
    subplot(1,2,1)
    errorbar(g_bins,KsnMod,ksnerrorinputs,'ko'); hold on
    xlabel('Unit ID'); ylabel('k_{sn}');
    SSz = Achi*KsnMod; %modeled long profile from ksn
    subplot(1,2,2)
    plot(S.distance,Sz,'.','color',[0.7 0.7 0.7]); hold on
    plot(S.distance,SSz,'k.');
end
%% Calculate Kmod erodibility values from Ksn and errors 
if isempty(p.Results.K)    
    En = ((E+E_error)-(E-E_error)).*rand(n,1)+(E-E_error);
    Ksn_n = zeros(n, length(KsnMod));
    Kn = zeros(n, length(KsnMod));
    % for loop through each ksn group to calculate n random values within
    % the error window
    for i = 1:length(KsnMod)
         t = ((KsnMod(i)+(KsnStd(i)))-(KsnMod(i)-KsnStd(i))).*rand(n,1)+(KsnMod(i)-KsnStd(i));
         Ksn_n(1:n,i) = t;
    end 
    % for loop though each iteration to calculate the K values from the
    % random Ksn and random E values
    for i = 1:n
        Kn(i,1:length(KsnMod)) = En(i)./Ksn_n(i,1:length(KsnMod));
    end
    Kn(1,:) = E./KsnMod; %the first iteration is with the best-fit K values
else
    if numel(K)>numel(KsnMod)
        error('There are more K values than geology types');
    end
    if numel(K)<numel(KsnMod)
        error('There are less K values than geology types');
    end
    Kn = zeros(n, length(K));
    % for loop through each k group to calculate n random values within
    % the error window
    for i = 1:length(K)
         t = ((K(i)+(K_error(i)))-(K(i)-K_error(i))).*rand(n,1)+(K(i)-K_error(i));
         Kn(1:n,i) = t;
    end 
    Kn(1,:) = K; %the first iteration is with the input K values
end
%% Run linear inversion block uplift variable K 
timesteps = NaN(50,n);
U_rates = NaN(50,n);
final_ts = zeros(n,1);
h = waitbar(0,strcat('Iterating inversion: 1/',string(n)));
figure(3); hold on
for i = 1:n
    Kmod_i = Kn(i,:);
    Kmod_i = Kmod_i(:);
    Stau = Achi*((Kmod_i).^-1);
    [U,tau_steps,Atau] = linear_inversion_block_uplift_variableK(DEM,Stau,tau_inc,'crita',crita,'mn',mn);
    tau_steps = tau_steps(:);
    timesteps(1:length(tau_steps),i) = tau_steps/1e6; 
    U_rates(1:length(U),i) = U*10^3;
    nts = rmmissing(timesteps(:,i)); %temporary list for final_ts
    final_ts(i) = nts(length(nts));
    stairs(rmmissing(timesteps(:,i)),rmmissing(U_rates(:,i)),'Color',[0.75 0.75 0.75]);
    waitbar(i/n,h,strcat('Iterating Inversion:',string(i+1),'/',string(n)));
end
close(h)
%calculate the average and median value at each timestep
Mean_Urate = zeros((max(max(timesteps)))./(tau_inc/1e6)+1,1);
Med_Urate = zeros((max(max(timesteps)))./(tau_inc/1e6)+1,1);
BF_timesteps = 0:(tau_inc/1e6):max(max(timesteps));
BF_timesteps = BF_timesteps(:);
for i = 1:length(Mean_Urate)
    Mean_Urate(i) = mean(rmmissing(U_rates(i,:)));
    Med_Urate(i) = median(rmmissing(U_rates(i,:)));
end
% calculate the median and mean final timestep to determine how long
% the best-fit record should extend
Med_TS = median(final_ts);
%if the median lies between two timesteps, round up to the nearest
if mod(Med_TS,(tau_inc/1e6)) ~= 0
    Med_TS = Med_TS + (mod(Med_TS,(tau_inc/1e6)));
end
% find the index of the median end and clip the best-fit record
ind = (Med_TS/(tau_inc/1e6))+1;
BF_timesteps = BF_timesteps(1:ind,:);
Med_Urate = Med_Urate(1:ind,:);
Mean_Urate = Mean_Urate(1:ind,:);
%Plot the data
figure(3)
stairs(BF_timesteps,Mean_Urate, 'Color', 'r', 'LineWidth',2);
stairs(BF_timesteps,Med_Urate, 'Color', 'b', 'LineWidth',2);
ylim([0 (max(Med_Urate)+0.25)]);
xlabel('time (Myr)'); ylabel('Uplift rate (mm/yr)');
annotation('textbox', [0.15,0.15,0.30,0.15],'String',['Red = Average U' newline 'Blue = Median' newline 'Gray = n iterations']);
hold off
%% Plot response time with input K values
if isempty(p.Results.K)
    K_tau = E./KsnMod;
    Stau = Achi*(K_tau.^-1);
    if isempty(p.Results.geol)
        figure(4); hold on
        scatter(Stau/1e6,Sz,2,[0.5 0.5 0.5])
        xlabel('Tau (Myr)'); ylabel('Elevation (m)');
        title('Response time uniform K'); hold off
    else
        figure(4); hold on
        scatter(Stau/1e6,Sz,2,S_discrete,'filled')
        colormap(parula);
        colorbar('Ticks',g_bins)
        xlabel('Tau (Myr)'); ylabel('Elevation (m)');
        title('Response time variable least-squares K'); hold off
    end
else
    Stau = Achi*((K).^-1);
    figure(4); hold on
    scatter(Stau/1e6,Sz,2,S_discrete,'filled')
    colormap(parula);
    colorbar('Ticks',g_bins)
    xlabel('Tau (Myr)'); ylabel('Elevation (m)');
    title('Response time variable K'); hold off
end
%% Plot Chi vs elevation
figure(5)
scatter(Schi,Sz,2,[0.5 0.5 0.5])
xlabel('\chi'); ylabel('Elevation (m)');
title('\chi plot'); hold off
end
% %%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions nessesary to complete the code above %%%%%
function [KsnMod,KsnStd,A,Schi] = invert_for_ksn(S,Sz,S_DA,S_discrete,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% invert_for_ksn_STREAMobj_prepped.m is a function based on data structures
% from TopoToolbox v. 2 (https://topotoolbox.wordpress.com/; Schwanghart
% and Scherler, 2014). Use of this function requires some familiarity with
% Matlab and TopoToolbox, but an example script is provided to help guide
% the user through several applications of the code.
%
% The code performs a linear inversion of the transformed stream
% channel distance variable chi (Perron and Royden, 2013) and channel
% elevation data over user-defined discrete intervals. The descrete
% intervals can be defined in a number of ways, for example, by a geologic
% map or by increments of stream channel elevation, distance, or chi. Some
% examples of how to discretize the channel network are given in the
% example script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% Date modified: 04/02/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%   Required [4]:
%       (1) S - STREAMobj from TopoToolbox
%       (2) Sz - Elevation of the stream network with respect to the outlet
%           (e.g. Sz = DEM.Z(S.IXgrid); Sz = Sz - min(Sz);)
%       (3) S_DA - Drainage area of the stream network in meters.
%       (4) S_discrete - spatial discretization of stream network. This
%           variable will be filled with a series of integers.
%
%   Optional [2]:
%       (1) mn - the m to n ratio or concavity used to calculate chi
%                {default = 0.5}
%       (2) inverse_option - technique used to invert the data.
%           Option 1: 'non_negative' - uses a non-negative least squares
%           inversion. This method provides 1-sigma uncertainties {default}
%           Option 2: 'Tikhonov' - uses a Tikhonov regularization to find 
%           in the least squares solution. Uncertainties cannot be provided
%           with this approach because the regulatization biases the model
%           to solve the inverse problem
%
% Outputs [4]:
%       (1) KsnMod - vector of ksn values for the least squares solution 
%                    ordered by integers in S_discrete.
%       (2) KsnStd - if the default non-negative least-squares inversion is
%                    used, this is a vector of 1-sigma uncertainties on the
%                    Ksn values calculated after accounting for serial
%                    correlation of residuals (e.g. Perron and Royden,
%                    2013). If the Tikhonov regularization is used, it is a
%                    vector of NaNs because this approach biases the
%                    problem precluding the ability to calculate
%                    uncertainties.
%       (3) A -      The forward model matrix used in the inversion. The
%                    forward matrix can be used to predict the theoretical
%                    steady-state elevation of the river work based on the
%                    KsnMod vector. (e.g. SSz = A*KsnMod).
%       (4) Schi -   the transformed distance variable, chi, for the
%                    STREAMobj
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Parse Inputs
p = inputParser;
p.FunctionName = 'invert_for_ksn';

% required inputs
addRequired(p,'S', @(x) isa(x,'STREAMobj')||isvector(x));
addRequired(p,'Sz', @(x) isvector(x));
addRequired(p,'S_DA',@(x) isvector(x));
addRequired(p,'S_discrete', @(x) isvector(x));

% optional inputs
addOptional(p,'mn', 0.5, @(x) isscalar(x));
addOptional(p,'inverse_option', []);

% declare variable names
parse(p, S, Sz, S_DA, S_discrete, varargin{:});
S   = p.Results.S;
Sz   = p.Results.Sz;
S_DA  = p.Results.S_DA;
S_discrete   = p.Results.S_discrete;
mn = p.Results.mn;

% some error handling
if length(S.distance) ~= length(S_discrete)
    error('The STREAMobj and S_discrete are note the same size')
elseif ~isempty(S_discrete(isnan(S_discrete)))
    error('S_discrete contains NaNs, which are not allowed');
elseif nanmin(Sz) ~= 0
    error('Sz is not calcaulted with respect to the outlet. Recalculate Sz as Sz = Sz - min(Sz)');
end


% inverse technique options
if isempty(p.Results.inverse_option)
    in_op = 1;
elseif strcmp(p.Results.inverse_option, 'non_negative')
    in_op = 1;
elseif strcmp(p.Results.inverse_option, 'Tikhonov')
    in_op = 2;
else
    error('fillOption is not "non_negative" or "Tikhonov"');
end

%% make the forward operator
[A,Schi] = MakeDiscreteChiAMatrix(S,S_DA,S_discrete,mn);

[N,q] =size(A);

% invert for ksn using preferred method.
if in_op == 1
    % invert for ksn.
    KsnMod = lsqnonneg(A,double(Sz));
    
    % compute the modeled elevations
    z_mod = A*KsnMod;
    
    % compute the residuals
    res = z_mod - Sz;
    
    % determine uncertainties on the mean ksn using a bootstrap approach
    % where the resampled size for each iteration is selected to minimize
    % serial correlation of residuals, a problem noted in Perron and Royden
    % (2013) regarding using the chi method to find ksn.
    %
    % The approach used here is bit ad hoc. It progressively thins the
    % dataset to define a percent data removed vs r-squared plot. This plot
    % generally follows a noisy trade curve. The code then fits a power law
    % function to the percent data removed vs r-squared data to find the
    % turning point in the data that is assumed to the be optimal amount of
    % thinning of the dataset to remove the impact of serial correlation.
    % The uncertainties are then determined using a bootstrap approach to
    % calculate ksn using the thinned the data and provide appropriate
    % 1-sigma uncertainties on the ksn estimates.
    
    % number of bootstrap runs
    n_runs = 5000;
    [KsnMod,KsnStd,~]=ksn_bootstrap(A,Sz,res,n_runs);
elseif in_op == 2
    % make q by q identity matrix
    I = eye(q,q);
    
    % calculate KsnPri
    KsnPri = ones(q,1)*(max(Sz)./max(Schi));
    
    % find the optimal dampening (e.g. smoothing) parameter (Gamma)
    Gam = logspace(-2,5,1000);
    MisFit = nan(size(Gam));
    
    for i = 1:length(Gam)
        KsnMod = KsnPri + (A'*A + Gam(i)^2*I)\A'*(Sz-A*KsnPri);
        MisFit(i) = (1/(N-q))*sqrt(sum((Sz - A*KsnMod).^2));
    end
    
    [ind,~] = turingPointFinder(1./Gam,MisFit);
    
    figure()
    plot(1./Gam,MisFit,'k-'); hold on
    plot(1/Gam(ind),MisFit(ind),'ko','markerfacecolor',[0.5 0.8 0.5]);
    xlabel('1/\Gamma'); ylabel('normalized misfit');
    H=text(1/Gam(ind)+0.002,MisFit(ind)+0.025,['best-fit \Gamma = ',num2str(Gam(ind))]);
    set(H,'Fontsize',10);
    title('Trade off curve')
    
    % now we can invert to find U following Goren eq (21) from Tarantola, 1987
    KsnMod = KsnPri + (A'*A + Gam(ind)^2*I)\A'*(Sz-A*KsnPri);
    
    % note, you can't get uncertainties on ksn with the Tikhonov
    % regularization because you bias the problem to find the least-squares
    % solution. In this case, we just fill the KsnStd vector with nans
    KsnStd = nan(size(KsnMod));
end



end
 %This function calculates chi and sets of the forward operator matrix
function [A,Schi] = MakeDiscreteChiAMatrix(S,S_DA,S_discrete,mn)
%
% MakeDiscreteChiAMatri.m will make an n by q matrix, A, that is n stream
% nodes long by q discrete classification units or domains. The matrix
% consists of the integral quantity chi traversed for each discrete domain
% such that sum(A,2) [sum of all the rows] will equal chi for the stream
% network. This matrix can be used to calculate average ksn per geologic
% map unit using a matrix inversion of elevation For example, ksnPerUnit =
% A\Sz, where Sz is river network elevation and ksnPerUnit is a vector that
% is q long with the average ksn of each discrete unit. Note that the above
% inversion allows for negative numbers, which isn't physically realiztic,
% so it is suggested to use the non-negative least-squares inverse command
% in Matlab, lsqnonneg, or some regulatization technique to invert the
% elevation and chi data for domain average ksn.
%
% Inputs:
% S             - TopoToolbox STREAMobj.
% DA            - Drainage area grid IN MAP UNITs (e.g. m^2) as a GRIDobj.
% S_discrete    - vector of STREAMobj nodes classified by use discretization.
% mn            - m/n or refence concavity (theta) value.
%
% Outputs:
% A             - A matrix as described above.
% Schi          - Chi for the stream network.
%
% Author: Sean F. Gallen
% Date modified: 06/08/2017
% email: sean.gallen{at}colostate,edu

p = inputParser;
p.FunctionName = 'MakeChiGeoAMatrix';
addRequired(p,'S', @(x) isa(x,'STREAMobj')||isvector(x));
addRequired(p,'S_DA', @(x) isvector(x));
addRequired(p,'S_discrete', @(x) isvector(x));
addRequired(p,'mn', @(x) isscalar(x));

parse(p,S,S_DA,S_discrete,mn);

% some error handling
if length(S.distance) ~= length(S_discrete)
    error('The STREAMobj and S_discrete are note the same size')
elseif ~isempty(S_discrete(isnan(S_discrete)))
    error('S_discrete contains NaNs, which are not allowed');
end


% get variables ready for chi integration
Schi = zeros(size(S.distance));
dSc = zeros(size(S.distance));
Six = S.ix;
Sixc = S.ixc;
Sx = S.distance;
Sa = S_DA.^-mn;

h = waitbar(0,'calculating \chi for full stream network...');
% calculating chi and tau_star for the entire river network
for lp = numel(Six):-1:1
    Schi(Six(lp)) = Schi(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
    dSc(Six(lp)) =  Schi(Six(lp))- Schi(Sixc(lp));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end
close(h);

% set up variables for the A matrix
discID = unique(S_discrete);

N = length(Schi);
q = length(discID);

A = zeros(N,q);

h = waitbar(0,'building the A matrix...');
for i = 1:q
    Dinds = find(S_discrete == discID(i));
    A(Dinds,i) = dSc(Dinds);
    for lp = numel(Six):-1:1
        A(Six(lp),i) = A(Sixc(lp),i) + A(Six(lp),i);
    end
    waitbar(i/q,h)
end
close(h)

end
function [ksn_mean,ksn_std,ksn_boot]=ksn_bootstrap(A,Sz,res,n_runs)
% ksn_bootstrap defines a resample size to minimize serial correlation of
% residuals and uses a bootstrap appoach to calculate the mean and standard
% deviation of ksn using a non-negative least squares regression.
%
% Inputs:
% (1) A:    matrix from "MakeDiscreteChiMatrix.m"
% (2) Sz:   stream network elevation matrix corresponding to A
% (3) res:  residual vector generated from non-negative least-sqaures
%           regression of A and Sz.
% (4) n_runs: number of iterations in the bootstrap. 1000 to 5000 is
%           probably good
%
% Outputs:
% (1) ksn_mean: mean ksn from bootstrap. This should be comparable to the
%               non-negative least-square regression result from A and Sz.
% (2) ksn_std:  standard devision of ksn from the bootstrap
% (3) ksn_boot: results from each bootstrap iteration.
%
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% date modified: 06/30/2019

p = inputParser;
p.FunctionName = 'ksn_bootstrap';
addRequired(p,'A', @(x) ismatrix(x));
addRequired(p,'Sz', @(x) isvector(x));
addRequired(p,'res', @(x) isvector(x));
addRequired(p,'n_runs', @(x) isscalar(x));

parse(p,A,Sz,res,n_runs);


%% vector of percentage of data removed
% This part is rather ad hoc, but seems to work well. Because the
% trade-off curve is caculated by bootstraping, the turing point varies by
% ~2-3% between model runs. If you don't find this statisfying, simply
% increase the number for 'iter' in the forloop below. I have chosen 200 to
% maximize the robustness of the "fit" and the computational time.
percent_removed = logspace(log10(1),log10(90),100)';
n_points = length(Sz) - round(length(Sz).*(percent_removed/100));
R = nan(size(percent_removed));
[Ar,Ac] = size(A);
ksn_trials = nan(Ac,n_runs);

% run a short bootstrap for each percentage to smooth out noisiness in the
% analysis
h = waitbar(0,'Finding optimal data thinning percentage...');
for i = 1:length(percent_removed)
    
    iter = 200;
    R_temp = nan(iter,1);
    for j = 1:iter
    % calculate the correlation chnages with thinning
    rand_inds = randi(Ar,n_points(i),1);
    Z_boot = Sz(rand_inds);
    A_boot = A(rand_inds,:);
    ksn_trials(:,i) =lsqnonneg(A_boot,double(Z_boot));
    
    SSz = A*ksn_trials(:,i);

    r_corr_temp = corrcoef(Sz, SSz);
    R_temp(j) = r_corr_temp(2);
    end
    R(i) = nanmean(R_temp);
    waitbar(i/length(percent_removed),h);
end
close(h)

%% smooth the data to make a better trade-off curve
inv_p_interp = linspace(min(1./percent_removed),max(1./percent_removed),1000)';
inv_R_interp = interp1(1./percent_removed,1./R,inv_p_interp,'pchip');
inv_R_interp = movavg(inv_R_interp,'simple',1000);

%% note, I will want to change this to a maximum curvature function
[idx,~] = turingPointFinder(inv_p_interp,inv_R_interp);

%% plot the trade-off curve
figure()
plot(1./percent_removed,1./R,'o','color',[0.7 0.7 0.7]); hold on
plot(inv_p_interp,inv_R_interp,'b-');
plot(inv_p_interp(idx),inv_R_interp(idx),'ko','markerfacecolor',[0.5 0.5 0.5]);
xlabel('1/Percentage of points removed'); ylabel('1/Correlation coefficient (R)')

H=text(inv_p_interp(idx)+(inv_p_interp(idx).*1e-3),inv_R_interp(idx)+...
    (inv_R_interp(idx).*5e-5),['percentage of points removed = ',...
    num2str(1./inv_p_interp(idx)),'%']);
set(H,'Fontsize',10);
title('Trade off curve')

% Use the "optimal" thinning percentage to bootstrap uncertainties.
% This will generally avoid the autocorrelation of residuals to give a
% better estimate of the true uncertainty.
rs_size = length(Sz) - round(length(Sz).*((1/inv_p_interp(idx))/100));

%% run the bootstrap
[Ar,Ac] = size(A);
ksn_boot = nan(Ac,n_runs);

h = waitbar(0,'Bootstrapping to get uncertainty...');

for i = 1:n_runs
    
    rand_inds = randi(Ar,rs_size,1);
    Z_boot = Sz(rand_inds);
    A_boot = A(rand_inds,:);
    ksn_boot(:,i) =lsqnonneg(A_boot,double(Z_boot));
    waitbar(i/n_runs,h);
end
close(h)
%% calculate mean and std of data.
ksn_mean = nanmean(ksn_boot,2);
ksn_std = nanstd(ksn_boot,0,2);
end
%% inversion function
function [Umod,tau_steps,A] = linear_inversion_block_uplift_variableK(DEM,Stau, tau_inc,varargin)
% Block uplift inversion code after Goren et al. (2014) JGR-Earth Surface.
% Inputs:
%       Required:
%       (1) DEM clipped to a watershed as a TopoToolbox GRIDobj
%       (2) Stau tau network calculated with K values for each geologic unit
%           assuming n = 1
%       (3) Tau increment in years.
%       (4) Stau - tau of the stream network 
%
%       Optional:
%       (1) crita - critical drainage area for channel head initiation in
%           m^2 (default = 1e6)
%       (2) mn - the m over n ratio (concavity) of river system (default =
%           0.5)
%       (3) flowOption: string of 'fill' or 'carve' for flow routing.
%          (optional) {default --> empty assumes that DEM is already 
%          filled or carved}
%       (4) Gam - Gamma dampening value 
%
% Outputs:
%       (1) A - forward model matrix
%       (2) Umod - recovered uplift history
%       (3) S - topotoolbox STREAMobj
%       (4) tau_steps - array of tau values at each step determined from
%       tau_inc
%     
% Author: Sean F. Gallen modified by James A. Fisher
% Date Modified: 05/18/2020
% email: sean.gallen[at]colostate.edu

%%Parse Inputs
p = inputParser;         
p.FunctionName = 'linear_inversion_block_uplift_variableK';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'Stau');
addRequired(p,'tau_inc');


% optional inputs
addOptional(p,'crita', 1e6, @(x) isscalar(x));
addOptional(p,'mn', 0.5, @(x) isscalar(x));
addOptional(p,'flowOption', []);
addOptional(p,'Gam',[],@(x) isscalar(x));

parse(p,DEM, Stau, tau_inc, varargin{:});
DEM   = p.Results.DEM;
Stau = p.Results.Stau;
tau_inc = p.Results.tau_inc;
crita = p.Results.crita;
mn = p.Results.mn;
Gam = p.Results.Gam;


%% Process grids
% set nan values if it hasn't already been done
DEM.Z(DEM.Z <= -9999) = NaN;

% declare cellsize
cs = DEM.cellsize;

% make direction FLOWobj
% flow routing options
if isempty(p.Results.flowOption)
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flowOption, 'fill')
    DEM = fillsinks(DEM);
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flowOption, 'carve')
    FD = FLOWobj(DEM,'preprocess','carve');
    DEM = imposemin(FD,DEM);
else
    error('fillOption is not "fill" or "carve"');
end

% make flow accumulation grid
%DA = flowacc(FD).*(FD.cellsize^2);

% Create stream network object (S)
S  = STREAMobj(FD,'minarea',crita/(cs)^2);
%S = trunk(S);% Include this command to run the inversion with only the
%trunk channel. Must also uncomment this command on line 181
S = klargestconncomps(S,1);

Sz = DEM.Z(S.IXgrid);



%% coding Goren et al., 2014 block uplift inversion model
% sort by elevation
table = sortrows([Sz,Stau]);
z_vec = table(:,1)-min(table(:,1)); %Creates base level of zero
tau_vec = table(:,2);

n = length(z_vec);  % number of nodes

% define desired time increments.
min_tau = floor(min(tau_vec));
max_tau = ceil(max(tau_vec));

if tau_inc > max_tau
    error(['Error: Your Tau increment is larger than the maximum Tau. ',...
        'The maximum Tau is: ', num2str(max_tau),' years. Reduce the the increment',...
        'so that you have ~5 to 10 tau increments between 0 and the maximum Tau.']);
end

% make time vector
%tau_steps = min_tau:tau_inc:max_tau; % time steps
tau_steps = min_tau : tau_inc : max_tau;
q = length(tau_steps);

% load the time matrix with data
A = zeros(n,q);
for i = 2:q
    inc_chi_t=tau_steps(i)-tau_steps(i-1);
    A(tau_vec >= tau_steps(i),i-1) = inc_chi_t;
end
for i = 1:n
    loc_sum = sum(A(i,:));
    loc = find(A(i,:)==0,1);
    A(i,loc) = tau_vec(i) - loc_sum; %this is the reminder
end

% calculate Upri
ZoverAsum = z_vec./sum(A,2);
ZoverAsum(ZoverAsum == Inf) = nan;
Upri = ones(q,1)*(1/n)*nansum(ZoverAsum);

% make q by q identity matrix
I = eye(q,q);

% impose dampening (e.g. smoothing) parameter, (Gamma)
if isempty(p.Results.Gam)
    Gam = logspace(-2,5,1000);

%% find the best-fit Gamma and plot results
    cols = jet(length(Gam));
    MisFit = nan(1,length(Gam));

    %figure()
    %subplot(2,2,1)
    for i = 1:length(Gam)
    % now we can invert to find U following Goren from Tarantola, 1987
        Umod = Upri + (A'*A + Gam(i)^2*I)\A'*(z_vec-A*Upri);
%     subplot(2,1,1);
%     stairs(chi_steps,Umod,'color',cols(i,:),'lineWidth',0.5); hold on
        MisFit(i) = (1/(n-q))*sqrt(sum((z_vec - A*Umod).^2));
    %subplot(2,1,2)
        %plot(1./Gam(i),MisFit(i),'.','color',cols(i,:)); hold on
    end

    [ind,~] = turingPointFinder(1./Gam,MisFit);
    %subplot(2,1,2)
    %plot(1/Gam(ind),MisFit(ind),'ko','markerfacecolor',[0.5 0.5 0.5]);
    %xlabel('1/\Gamma'); ylabel('normalized misfit');

    % define best-fit Gamma and do the final inversion
    Gam = Gam(ind);
else
end
%% invert for uplift history with best damping factor
Umod = Upri + (A'*A + Gam^2*I)\A'*(z_vec-A*Upri);

% plot uplift results
%subplot(2,2,2)
% incremental uplift
%stairs(tau_steps./1e6,Umod.*1000,'color',[0 0 0],'lineWidth',2); hold on
%xlabel('\tau (Myr)'); ylabel('U (mm/yr)');

% cummulative uplift
%upWRTtime = cumtrapz(tau_steps,Umod);
%subplot(2,2,4)
%stairs(tau_steps./1e6,upWRTtime,'color',[0 0 0],'lineWidth',1); hold on
%xlabel('\tau (Myr)'); ylabel('total baselevel fall(m)');

% plot observed and best-fit model results
%subplot(2,2,3)
%plot(tau_vec./1e6,z_vec,'.','color',[0.5 0.5 0.5]); hold on
%plot(tau_vec./1e6,A*Umod,'k.');%,'color',Pcols(k,:));
%xlabel('\tau (Myr)'); ylabel('elevation (m)')
%legend({'observed','modeled'});
end

function [idx,perpDist] = turingPointFinder(x,y)
% function finds "corner" on a plot with x and y asymptotes
%
% Author: Sean F. Gallen
% Date modified: 03/23/2016
% email: sean.gallen{at}colostate.edu

% Two endpoints on the curve "data"
xend = [x(1) x(end)];
yend = [y(1) y(end)];
% The slope of the line connecting the two endpoints
m = ( yend(2) - yend(1) )/( xend(2) - xend(1) );
% Point on the curve (xc,yc), point on the line (xl,yl)
perpDist = zeros(length(x),1);
for i = 1:length(x)
    xc = x(i) ; yc = y(i);
    yl = ( (m * xc) + (m^2 * yc) - (m * xend(1)) + yend(1) )/(1+ m^2);
    xl = xc - m*(yl - yc);
    % distance^2
    d2 = (xl - xc)^2 + (yl - yc)^2;
    perpDist(i) = (d2);
end
[~, idx] = max(perpDist);
end
