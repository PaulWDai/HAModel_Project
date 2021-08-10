function [K_tr,r_tr] = TransitionPath(tau,tau_tr,time_tr,r,w,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: 
%   Solve the transition path after a policy shock (capital income tax
%   rate from 0 to 20%, for example)

% Method:
%   Solve out the initial and ending stationary equilibrium
%   Guess a length of transition and a path of aggregate capital level
%   Simulate the transition (density of different state agents)
%   Check the convergence of aggregate path and make sure the density is
%   similar in the last two periods. Otherwise, keep iteration or adjust
%   the length of transition.

% Input:
%   tau: capital income tax, in  this quesiton, it is set to 0
%   tau_tr: capital income tax during transition
%   time_tr: length of transition
%   r,w: initial guess of r,w for solving stationary equilibrium
%   p:  stationary distribution of markov state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global alpha beta sigma delta ...
    infty ygrid kgrid Kgrid sgrid nygrid nkgrid nsgrid nKgrid ...
    iter_max err_max tm_y  rho sig_y...
    Kmin Kmax ...
    obs time_win pre_burn ...
    fig_location options


% steady state: pre-reform
rw1 = ...
    fsolve(@(rw_vec)solve_rw_markov(rw_vec,rho,sigma,sig_y,tau),[r,w],options);
[~,k_density1,agg_k1,opt_k_idx1,density1] = ...
    solve_rw_markov(rw1,rho,sigma,sig_y,tau);

% steady state: post-reform
rw2 = ...
    fsolve(@(rw_vec)solve_rw_markov(rw_vec,rho,sigma,sig_y,tau_tr),[r,w],options);
[~,k_density2,agg_k2,opt_k_idx2,density2] = ...
    solve_rw_markov(rw2,rho,sigma,sig_y,tau_tr);

% initial guess
err_tr = 10^6;
iter = 0;
iter_tr = 300;
err_tr_max = 10^-2;
err_density_max = 10^-2;
time_tr_vec = 1:1:time_tr;
K_tr = (agg_k2 - agg_k1)/(time_tr-1)*(time_tr_vec-1)+agg_k1;

while iter<iter_tr && err_tr > err_tr_max
iter = iter+1;

L_tr = sum(p.*ygrid);
w_tr = (1-alpha)*L_tr^(-alpha)*K_tr.^(alpha);
r_tr = alpha*L_tr^(1-alpha)*K_tr.^(alpha-1)-delta;
indT_tr = tau*r_tr.*K_tr/1;

opt_k_tr = nan(nygrid,nkgrid,time_tr);
tm2D_tr = nan(nygrid*nkgrid,nygrid*nkgrid,time_tr);
density_tr = nan(nygrid*nkgrid,time_tr);

for time = 1:time_tr
    r = r_tr(time);
    w = w_tr(time);
    tau = .2;
    ind_T = indT_tr(time);
    [opt_k_idx,~] = ValueFunInfinite(r,w,rho,sigma,sig_y,tau,ind_T);
    opt_k_tr(:,:,time) = opt_k_idx;
    
    % transitional matrix of (y,k) pair
    sim_tm2D = TranMatrixEconomy(opt_k_idx,tm_y);
    tm2D_tr(:,:,time) = sim_tm2D;
    
    if time == 1
        density_last = reshape(density1,[nygrid*nkgrid,1]);
    else
        density_last = density_tr(:,time-1);
    end
        density_last = density_last';
        density_temp = density_last*sim_tm2D;
        density_tr(:,time) = density_temp';
end
density_tr3D = reshape(density_tr,[nygrid,nkgrid,time_tr]);
K_tr_new = permute(sum(density_tr3D,1),[2,3,1]);
K_tr_new = kgrid*K_tr_new;

err_tr = max(abs(K_tr_new-K_tr),[],'all')
K_tr = 0.3*K_tr_new + 0.7*K_tr ;
end

err_density = max(abs(density2-density_tr(:,time_tr)));
if err_density > err_density_max
    disp('Extend the Period of Transition')
else
    if err_tr < err_tr_max
        disp('Transition Path is Solved')
    else
        disp('No Covergence on Aggregate Capital')
    end
end

r_tr = alpha*K_tr.^(alpha-1)*L_tr^(1-alpha)-delta;

figure();
%subplot(2,1,1);
plot(1:time_tr,K_tr,'LineWidth',3.0);
xlabel('Time'),ylabel('$K$');
title('Transition Path of $K$');

% remark: the after-depreciation rental price of capital lies in
% (-delta,1/beta-1)

%{
subplot(2,1,2);
plot(1:time_tr,r_tr,'LineWidth',3.0);hold on;
plot(1:time_tr,-delta*ones(1,time_tr),'LineWidth',2.0);hold on;
plot(1:time_tr,(1/beta-1)*ones(1,time_tr),'LineWidth',2.0);
xlabel('Time'),ylabel('$r$');
legend('$r$','$-\delta$','$\rho$');title('Transition Path of $r$');
%}
saveas(gcf,fig_location+'Fig_Sec3_2.png')


end

