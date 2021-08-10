function [mean_c,mean_k,std_c,std_k] = ...
    solve_PartialLifeCycle(beta0,sigma0,rho0,sig_y0,inc,con,survs0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Given a set of parameters beta, sigma, sig_y, solve the partial life
%   cycle model.

% Output:
%   Life cycle profile: mean_c, mean_k, std_c, std_k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global w r beta sigma rho sig_y theta ...
    nygrid nwk nrt nkgrid nagrid ...
    awkgrid artgrid agrid kgrid amin amax abar ss ...
    iter_max err_max mu infty

% in case of ambiguity using global variables

beta = beta0;
sigma = sigma0;
rho = rho0;
sig_y = sig_y0;
survs = survs0;
% partial life cycle parameters
w = 1;              % wage
r = .01;            % return to capital
% beta = 1/(1+.04);   % discount
% sigma = 3;          % CRRA utility
% rho = .9;           % income process autocorrelation
mu = 0;             % mean of random shock 
% sig_y = .4;         % income process variance
theta = .5;         % ratio of social security income to income in age 45
infty = 999999;


% grid setting
nygrid = 7;

nkgrid = 250;       % grid for asset
kmin = 0;
kmax = 200;
kgrid = linspace(kmin,kmax,nkgrid);

amin = 21;          % grid for age
abar = 65;          % the last age of work
amax = 80;

nagrid = amax-amin+1;
nwk = abar-amin+1;      % number of age of work
nrt = amax-(abar+1)+1;  % number of age of retired
agrid = linspace(amin,amax,nagrid);

awkgrid = amin:1:abar;      % age of work grid
artgrid = abar+1:1:amax;    % age of retire grid

ss = theta*inc(length(inc),1);  % social security

% iteration setting
iter_max = 1000;
err_max = 10^-6;

% simulation
rng('default')
obs = 10000;   % sample size

%% Partial Equilibrium Life Cycle Model
% Tauchen
[tm_wk_y,ygrid] = ...
     Tauchen(nygrid,rho,sig_y);

% Transitional Matrix
[tm_wk,tm_rt] = ...
    TransMatrixPartial(tm_wk_y);

% Period utility
[mutilrt,mutilwk,Vterminal] = ...
    PeriodUtilityPartialLifeCycle(ygrid,inc);

% Value function iteration
[opt_krt_idx,opt_kwk_idx,opt_cwk,opt_crt] = ...
        ValFunIterPartialLifeCycle(tm_rt,tm_wk,mutilrt,mutilwk,Vterminal,inc,ygrid,survs);
%{
figure();
subplot(2,1,1);
s = surf(opt_kwk_idx);view(0,90);s.EdgeColor = 'none';
subplot(2,1,2);
s = surf(opt_krt_idx);view(0,90);s.EdgeColor = 'none';
%}
        
% Simulation

[mean_c,mean_k,std_c,std_k] = ...
    Simulation(obs, opt_krt_idx, opt_kwk_idx, opt_crt, opt_cwk, ygrid, inc,tm_wk_y);


end

