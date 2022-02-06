%% Get ksn distribution with errors
% this script is will run the "invert_for_ksn.m" function to find the
% least-squares ksn for each rock-type and output a table of the ksn
% distrubution.

% dist = [Rock-type ID; ksn; ksn error; # datapoints (stream nodes)]

% This script is meant to be used for comparisons of uplift histories
% scross multiple catchments. The same K values must be used for the same
% rock type. The ksn distributions can be used to calculate a weighted
% average K values for each rock-type across multiple catchments.

%Author:James A. Fisher jamesfisher.gis@gmail.com
%Sub-functions by: Sean F. Gallen; sean.gallen@colostate.edu
%
%Last Modified: 4/27/21
%% load data
mn = 0.45; %SET THETA. Theta must be consistent across all catchments
DEM = GRIDobj('name_ws.txt');
geol = GRIDobj('name_geol.txt');
%% Prep data
FD = FLOWobj(DEM);
DA = flowacc(FD);
DA = DA.*(FD.cellsize^2);
S  = STREAMobj(FD,'minarea',1e6/(FD.cellsize^2));
S = klargestconncomps(S,1);
Sz = DEM.Z(S.IXgrid);
Sz = Sz - min(Sz);
S_DA = DA.Z(S.IXgrid);
%% Binning by geology
S_discrete = geol.Z(S.IXgrid); %geology for each point on the stream network.
g_bins = unique(S_discrete); %Bin values set to unique values in S_discrete
% run the non-negative least-squares inversion and plot results
[KsnMod,KsnStd,A,Schi] = invert_for_ksn(S,Sz,S_DA,S_discrete,'mn',mn);
% plot the data
figure(103)
title('least-square ksn by geologic unit')
subplot(1,2,1)
errorbar(g_bins,KsnMod,KsnStd,'ko'); hold on
xlabel('Unit ID'); ylabel('k_{sn}');
SSz = A*KsnMod;
subplot(1,2,2)
plot(S.distance,Sz,'.','color',[0.7 0.7 0.7]); hold on
plot(S.distance,SSz,'k.');
%% Calculate distribution of ksn values. Output gets written to dist. 
[GC,GR] = groupcounts(S_discrete); 
dist = [GR,KsnMod,KsnStd,GC]; %distribution of ksn values [Rock-type ID; ksn; ksn error; # datapoints (stream nodes)]
%% Plot long profile colored with geologic unit 
figure(1)
scatter(S.distance,Sz,2,S_discrete,'filled')
colorbar('Ticks',g_bins)
title('long profile colored by geologic unit')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copied functions to complete the commands above
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
% Example:
% %% load and prep the data
% DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
% DEM = fillsinks(DEM);
% FD = FLOWobj(DEM);
% DA = flowacc(FD).*(FD.cellsize^2);
% S  = STREAMobj(FD,'minarea',1e6/(FD.cellsize^2));
% S = klargestconncomps(S,1);
% Sz = DEM.Z(S.IXgrid);
% Sz = Sz - min(Sz);
% S_DA = DA.Z(S.IXgrid);
% 
% %% Bin the data by elevation increments
% % set bin size arbitrarily to 200 m
% bin_size = 200;
% min_z = floor(nanmin(Sz)/bin_size)*bin_size;
% max_z = ceil((nanmax(Sz))/bin_size)*bin_size;
% 
% z_bins1 = min_z:bin_size:max_z;
% 
% bin_ids = 1:length(z_bins1)-1;
% z_bins = nan(size(bin_ids));
% 
% S_discrete = zeros(size(Sz));
% 
% for i = 1:length(bin_ids)
%     S_discrete(Sz >= z_bins1(i) & Sz < z_bins1(i+1)) = bin_ids(i);
%     z_bins(i) = (z_bins1(i)+z_bins1(i+1))/2;
% end
% 
% %% run the non-negative least-squares inversion and plot results
% [KsnMod,KsnStd,A,Schi] = invert_for_ksn(S,Sz,S_DA,S_discrete);
% % plot the data
% figure(99)
% subplot(1,2,1)
% errorbar(z_bins,KsnMod,KsnStd,'ko'); hold on
% SSz = A*KsnMod;
% subplot(1,2,2)
% plot(S.distance,Sz,'.','color',[0.7 0.7 0.7]); hold on
% plot(S.distance,SSz,'k.');
% 
% %% run inversion with Tikhonov regularization and plot results
% [KsnMod,KsnStd,A,Schi] = invert_for_ksn(S,Sz,S_DA,S_discrete,'inverse_option','Tikhonov');
% figure(99)
% subplot(1,2,1)
% plot(z_bins,KsnMod,'bo');
% xlabel('Mean streamwise distance (m)'); ylabel('k_{sn}');
% SSz = A*KsnMod;
% subplot(1,2,2)
% plot(S.distance,SSz,'b.');
% xlabel('Distance (m)'); ylabel('Elevation (m)');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%% additional functions used in code above %%%%%%%%%%%%%%%%%%%%

%% This function calculates chi and sets of the forward operator matrix
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

%% This function finds the turning point of the misfit function
function [idx,perpDist] = turingPointFinder(x,y)
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% date modified: 05/03/2016

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
