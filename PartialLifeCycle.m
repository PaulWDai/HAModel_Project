% Project 2
% Weifeng Dai

clear;
clc;
close;
%% Parameterization
global w r beta sigma rho sig_y theta ...
    nygrid nwk nrt nkgrid nagrid ...
    awkgrid artgrid agrid kgrid amin amax abar ss ...
    iter_max err_max mu infty
% partial life cycle parameters
w = 1;              % wage
r = .01;            % return to capital
beta = 1/(1+.04);   % discount
sigma = 3;          % CRRA utility
rho = .9;           % income process autocorrelation
mu = 0;             % mean of random shock 
sig_y = .4;         % income process variance
theta = .5;         % ratio of social security income to income in age 45
infty = 999999;
% read files of life cycle data;
incfile = fopen('incprofile.txt','r');
consfile = fopen('consprofile.txt','r');
survsfile = fopen('survs.txt');
formatSpec = '%f';
inc = fscanf(incfile,formatSpec);
con = fscanf(consfile,formatSpec);
con = reshape(con,[2,.5*length(con)])';

for i = 22:87
    cons(1,i-21) = mean(con(4*(i-22)+1:4*(i-22)+4,2));
end

survs = fscanf(survsfile,formatSpec);
survs = survs(length(survs)-nagrid+1:length(survs),1);

% grid setting
nygrid = 7;

nkgrid = 200;       % grid for asset
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
obs = 5000;   % sample size

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





