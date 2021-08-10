% General.m
% Project 2: Section 3 General Equilibrium

clear;
clc;
close;

%% Parameterization
global alpha beta sigma delta ...
    infty ygrid kgrid Kgrid sgrid nygrid nkgrid nsgrid nKgrid ...
    iter_max err_max tm_y ab bb ag bg rho sig_y...
    Kmin Kmax ...
    obs time_win pre_burn ...
    fig_location options K_ss

% partial life cycle parameters

beta = 1/(1+.04);   % discount factor
sigma = 3;          % CRRA utility
rho = .9;           % income process autocorrelation
mu = 0;             % mean of random shock
sig_y = .4;         % income process variance
alpha = .36;        % capital share
delta = .08;        % depreciation
infty = 99999;
% grid setting
nygrid = 7;

nkgrid = 200;       % grid for asset
kmin = 0.1;
kmax = 80;
kgrid = linspace(log(kmin),log(kmax),nkgrid); % individual capital
kgrid = exp(kgrid);

tm_s = [.9,.1;.1,.9];
sgrid = [0.95,1.05];
nsgrid = length(sgrid);

obs = 1000;
time_win = 200;
pre_burn = round(1/4*time_win);

nKgrid = 100;   % average capital per efficient labor
% steady state w/o uncertainty value
Kmin = 0.35*log(kmax);
Kmax = .5*log(kmax);
Kgrid = linspace(Kmin,Kmax,nKgrid);

% iteration setting
iter_max = 1000;
err_max = 10^-5;

% simulation
rng('default');

% graph setting
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
fig_location = "/Users/weifengdai/Dropbox/My Mac (MacBook Pro)/Documents/Undergrad/2021Spring/Heterogeneous Macro/Project/Project2/figure/";
options = optimset('Algorithm', 'levenberg-marquardt','Display','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 1. Stationary Distribution %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 (a) Benchmark: Wealth and Consumption Inequality

[tm_y,ygrid] = Tauchen(nygrid,rho,sig_y);
p = StationaryStatePDF(tm_y);

% initialization

r = 0.01;w = 1;

tau = 0;    % capital income tax
ind_T = 0;  % individual level transfer from gov.


disp('-----------------------------------------');
disp('Stationary Distribution Using Markov')
tic
rw_markov = ...
    fsolve(@(rw_vec)solve_rw_markov(rw_vec,rho,sigma,sig_y,tau),[r,w],options);
[~,k_density_markov,agg_k_markov,~,~,sample_c,sample_k] = ...
    solve_rw_markov(rw_markov,rho,sigma,sig_y,tau);
toc

disp('-----------------------------------------');
disp('Stationary Distribution Using Monte Carlo');
tic
rw_monte_carlo= ...
    fsolve(@(rw_vec)solve_rw_monte_carlo(rw_vec,rho,sigma,sig_y,tau,ind_T),[r,w],options);
[~,k_density_monte_carlo,agg_k_monte_carlo] =  ...
    solve_rw_monte_carlo(rw_monte_carlo,rho,sigma,sig_y,tau,ind_T);
toc


% Plot the result
Fig_Sec3_1_a(k_density_markov,k_density_monte_carlo,sample_k,sample_c);


%% 1 (b) Parameter Investigation


% use markov considering the time use and accuracy
vsigma = [1.01,3,5];
vrho = [0.3,0.5,0.7];
vsig_y = [.3,.5,.7];

[r_vsigma,w_vsigma,agg_k_vsigma,k_density_vsigma,...
          r_vrho,w_vrho,agg_k_vrho,k_density_vrho,...
          r_vsig_y,w_vsig_y,agg_k_vsig_y,k_density_vsig_y] = ...
          StationaryParameters(vsigma,vrho,vsig_y,tau);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% 2. Transition Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_tr = 70;   % period for transition
tau_tr = .2;    % policy shock: capital income tax

tic
[K_tr,r_tr] = TransitionPath(tau,tau_tr,time_tr,r,w,p);
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 3. With Aggregate Risk %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (a) Compute the approximated equilibrium with boundedly rational agents

%{
% law of motion
% baseline from Krusell and Smith (1998 JPE)

ag0 = 0.100;bg0 = 0.961;ab0 = 0.095;bb0 = 0.961;
coef = [ag0,bg0,ab0,bb0];

iter_KS = 0;err_KS = 10^6;

while iter_KS < iter_max && err_KS > 10^-3
    iter_KS = iter_KS+1;
    % solve a Krusel and Smith model given (a,b). 
    [coef_new,ols_g_Rsquared,olg_b_Rsquared,business_stat] = ...
        KruselSmith(coef,tm_y,tm_s,p);
    err_KS = max(abs(coef_new-coef));
    coef = coef_new
end

disp(coef);
disp(business_stat);
disp(ols_g_Rsquared,olg_b_Rsquared);

%}
