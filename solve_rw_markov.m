function [output,sim_k_density,agg_k_new,opt_k_idx,sim_density,sample_c,sample_k] =...
    solve_rw_markov(rw_vec,rho,sigma,sig_y,tau)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: 
%   Find the stationary general equilibrium of infinitely lived heterogeneous
%   agent economy with labor productivity shock

% Method:
%   Given wage rate and rental price of capital, solve out the policy
%   functions and get the transitional matrix of entire economy.
%   Solve the stationary distribution of the state of entire economy, and
%   Get the measure of agents with different states.
%   Update the wage rate and rental price of capital until convergence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = rw_vec(1);
w = rw_vec(2);

global alpha beta delta ...
    infty ygrid kgrid nygrid nkgrid ...
    iter_max err_max

[tm_y,ygrid] = Tauchen(nygrid,rho,sig_y);

% Markov
p_final = StationaryStatePDF(tm_y);
agg_l = sum(p_final.*ygrid);
agg_k = (alpha/(r+delta))^(1/(1-alpha))*agg_l;

ind_T = tau*r*agg_k;

[opt_k_idx,opt_c,~] = ValueFunInfinite(r,w,rho,sigma,sig_y,tau,ind_T);

sim_tm2D = TranMatrixEconomy(opt_k_idx,tm_y);

iter_sim_markov =0;
err_sim_markov = 10^6;
% initial guess
sim_density = ones(1,length(sim_tm2D))*1/length(sim_tm2D);

while iter_sim_markov < iter_max && err_sim_markov>err_max
    iter_sim_markov = iter_sim_markov+1;
    sim_density_new = sim_density*sim_tm2D;
    err_sim_markov = max(abs(sim_density_new-sim_density));
    sim_density = sim_density_new;
end



sim_density = reshape(sim_density,[nygrid,nkgrid]);
sim_k_density = sum(sim_density,1);

agg_k_new = sum(sim_k_density.*kgrid);
agg_l_new = sum(p_final.*ygrid);

rnew = alpha*agg_k_new^(alpha-1)*agg_l_new^(1-alpha)-delta;
wnew = (1-alpha)*agg_l_new^(-alpha)*agg_k_new^(alpha);

output(1) = abs(rnew-r);
output(2) = abs(wnew-w);

opt_c_vec = reshape(opt_c,[nygrid*nkgrid,1]);
opt_k_vec = reshape(kgrid(opt_k_idx),[nygrid*nkgrid,1]);
sim_density = reshape(sim_density,[nygrid*nkgrid,1]);

repeat = round(sim_density*10000);
count = 1;
for i = 1:length(opt_c_vec)
    count_end = count+repeat(i,1);
    sample_c(count:count_end-1,1) = opt_c_vec(i,1);
    sample_k(count:count_end-1,1) = opt_k_vec(i,1);
    count = count + repeat(i,1);
end

end

