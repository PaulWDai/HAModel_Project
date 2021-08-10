function [mean_c,mean_k,std_c,std_k] = ...
    Simulation(obs, opt_krt_idx, opt_kwk_idx, opt_crt, opt_cwk, ygrid,inc,tm_wk_y)

global w r rho sig_y  ...
    nygrid nwk nrt nkgrid nagrid...
    agrid kgrid amin  ss ...
    mu

rng('default');
%% Vectorization
% to use sub2ind to reduce the hierachy of loops
vopt_krt_idx = reshape(opt_krt_idx,[1,(nrt-1)*nkgrid]);
vopt_kwk_idx = reshape(opt_kwk_idx,[1,nwk*nygrid*nkgrid]);
vopt_crt = reshape(opt_crt,[1,(nrt-1)*nkgrid]);
vopt_cwk = reshape(opt_cwk,[1,nwk*nygrid*nkgrid]);

sim_idx = nan(obs,nagrid);
% sim_epsilon = nan(obs,nagrid);
sim_k = nan(obs,nagrid);
sim_c = nan(obs,nagrid);
sim_a = repmat(agrid,[obs,1]);
sim_a_idx = sim_a - amin + 1;

%% Simulation Initialization

tm_y_markov = dtmc(tm_wk_y);
sim_y_idx = simulate(tm_y_markov,nwk,'X0',[0,0,0,obs,0,0,0])';
sim_y = ygrid(sim_y_idx);
sim_k(:,1) = 0; % a0 = 0;

%{

% Here is the wrong version I wrote before!!!

% Firstly, I did not fully utilize the discretized y; But it is not
techanically a mistake;

sim_epsilon = normrnd(mu,sig_y,size(sim_epsilon));

for period = 2:nagrid
    sim_y(:,period) = rho*sim_y(:,period-1) ...
        + sim_epsilon(:,period-1);
end

% Secondly, the simulated y index using the methodology below is not
correct, which gives a wrong version of life cycle performance.
% The ygrid is not equally spaced. I should use find to determine the
related grid index.

sim_y_idx = round((sim_y-min(ygrid))/(max(ygrid)-min(ygrid))*(nygrid-1))+1;
sim_y_idx(sim_y_idx<1) = 1;
sim_y_idx(sim_y_idx>nygrid) = nygrid;
%}

sim_k_idx = sim_k;
sim_k_idx(:,1) = round((sim_k(:,1)-min(kgrid))/(max(kgrid)-min(kgrid))*(nkgrid-1))+1;

for age = 1:nwk-1
    idx = sub2ind([nygrid,nwk,nkgrid],sim_y_idx(:,age),sim_a_idx(:,age),sim_k_idx(:,age));
    sim_idx(:,age) = idx;
    sim_k_idx(:,age+1) = vopt_kwk_idx(idx);
    sim_c(:,age+1) = vopt_cwk(idx);
    sim_k_idx(sim_k_idx>nkgrid) = nkgrid;
    sim_k_idx(sim_k_idx<1) = 1;
end

for age = nwk:nwk+nrt-2
    idx = sub2ind([nrt-1,nkgrid],sim_a_idx(:,age)-nwk+1,sim_k_idx(:,age));
    sim_idx(:,age) = idx;
    sim_k_idx(:,age+1) = vopt_krt_idx(idx);
    sim_c(:,age+1) = vopt_crt(idx);
end

sim_k_idx(:,nwk+nrt) = 1; % no asset at all
sim_k = kgrid(sim_k_idx);

sim_c(:,1) = w*sim_y(:,1).*inc(1) + (1+r)*sim_k(:,1)-sim_k(:,2);
sim_c(:,nagrid) =  (1+r)*sim_k(:,nagrid) + ss;

mean_c = mean(sim_c,1);
std_c = std(sim_c,1);
mean_k = mean(sim_k,1);
std_k = std(sim_k,1);

end

