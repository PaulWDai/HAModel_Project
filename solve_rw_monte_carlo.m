function [output,sim_k_density_monte_carlo,agg_k] = ...
    solve_rw_monte_carlo(rw_vec,rho,sigma,sig_y,tau,ind_T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: 
%   Find the stationary general equilibrium of infinitely lived heterogeneous
%   agent economy with labor productivity shock

% Method:
%   Given wage rate and rental price of capital, solve out the policy
%   functions and simulate a sample of households
%   Get the equilibrium after a long time window, where the economy
%   converges to the stationary distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('default');

global alpha beta delta ...
    infty ygrid kgrid  nygrid nkgrid  ...
    tm_y obs time_win

r = rw_vec(1);
w = rw_vec(2);

[tm_y,ygrid] = Tauchen(nygrid,rho,sig_y);


markov = dtmc(tm_y);
sim_y_idx = nan(obs,time_win);

for i = 1:obs
    v1 = simulate(markov,time_win-1)';
    sim_y_idx(i,:) = v1;
end


[opt_k_idx,~] = ValueFunInfinite(r,w,rho,sigma,sig_y,tau,ind_T);
vopt_k_idx = reshape(opt_k_idx,[nygrid*nkgrid,1]);

% Simulation
sim_k_idx = nan(obs,time_win);
sim_k_idx(:,1) = ones(obs,1);

for j = 1:time_win-1
   ind = sub2ind([nygrid,nkgrid],sim_y_idx(:,j),sim_k_idx(:,j));
   sim_k_idx(:,j+1) = vopt_k_idx(ind);
end

sim_k = kgrid(sim_k_idx);
agg_k = sum(sim_k,1);
sim_y = ygrid(sim_y_idx);
agg_l = sum(sim_y,1);

agg_k = agg_k(time_win);
agg_l = agg_l(time_win);
rnew = alpha*agg_k^(alpha-1)*agg_l^(1-alpha)-delta;
wnew = (1-alpha)*agg_l^(-alpha)*agg_k^alpha;

% plot the distribution of 
[sim_k_density,~] = histcounts(sim_k(:,time_win),kgrid);
sim_k_density = sim_k_density/sum(sim_k_density);
sim_k_density_monte_carlo = [1-sum(sim_k_density),sim_k_density];

output(1) = abs(rnew-r);
output(2) = abs(wnew-w);

end

