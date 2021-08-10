function [r_vsigma,w_vsigma,agg_k_vsigma,k_density_vsigma,...
          r_vrho,w_vrho,agg_k_vrho,k_density_vrho,...
          r_vsig_y,w_vsig_y,agg_k_vsig_y,k_density_vsig_y] = ...
          StationaryParameters(vsigma,vrho,vsig_y,tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: 
%   Comparison of equilibrium results with different sets of parameters
%   agent economy with labor productivity shock

% Method:
%   For a fixed group of parameter
%   Given wage rate and rental price of capital, solve out the policy
%   functions and get the transitional matrix of entire economy.
%   Solve the stationary distribution of the state of entire economy
%   Get the measure of agents with different states.
%   Plot the result

% Input:
%   vector of parameter candidate sigma, rho, sig_y, where v denotes vector
%   tau, capital income tax, in  this quesiton, it is set to 0

% Output:
%   r, w, agg_k,k_density: rental price of capital, wage, aggregate
%   capital, distribution of capital under different pairs of parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global rho sigma sig_y fig_location options
r = .01;
w = 1;

for i = 1:length(vsigma)
    rw_markov = ...
    fsolve(@(rw_vec)solve_rw_markov(rw_vec,rho,vsigma(i),sig_y,tau),[r,w],options);

    [~,k_density_markov,agg_k_markov,~,~,~] = ...
    solve_rw_markov(rw_markov,rho,vsigma(i),sig_y,tau);
    r_vsigma(i) = rw_markov(1); 
    w_vsigma(i) = rw_markov(2);
    agg_k_vsigma(i) = agg_k_markov;
    k_density_vsigma(i,:) = k_density_markov;
end


for i = 1:length(vrho)
    rw_markov = ...
    fsolve(@(rw_vec)solve_rw_markov(rw_vec,vrho(i),sigma,sig_y,tau),[r,w],options);
    [~,k_density_markov,agg_k_markov,~,~,~] = ...
    solve_rw_markov(rw_markov,vrho(i),sigma,sig_y,tau);
    r_vrho(i) = rw_markov(1); 
    w_vrho(i) = rw_markov(2);
    agg_k_vrho(i) = agg_k_markov;
    k_density_vrho(i,:) = k_density_markov;
end


for i = 1:length(vsig_y)
    rw_markov = ...
    fsolve(@(rw_vec)solve_rw_markov(rw_vec,rho,sigma,vsig_y(i),tau),[r,w],options);
    [~,k_density_markov,agg_k_markov,~,~,~] = ...
    solve_rw_markov(rw_markov,rho,sigma,vsig_y(i),tau);
    r_vsig_y(i) = rw_markov(1); 
    w_vsig_y(i) = rw_markov(2);
    agg_k_vsig_y(i) = agg_k_markov;
    k_density_vsig_y(i,:) = k_density_markov;
end


% plot the result


figure();
subplot(3,3,1);
plot(vsigma,r_vsigma,'r','LineWidth',2.0);
xlabel('$\sigma$');ylabel('$r$')
subplot(3,3,2);
plot(vsigma,w_vsigma,'b','LineWidth',2.0);
xlabel('$\sigma$');ylabel('$w$')
subplot(3,3,3);
plot(vsigma,agg_k_vsigma,'g','LineWidth',2.0);
xlabel('$\sigma$');ylabel('$K$')

subplot(3,3,4);
plot(vrho,r_vrho,'r','LineWidth',2.0);
xlabel('$\rho$');ylabel('$r$')
subplot(3,3,5);
plot(vrho,w_vrho,'b','LineWidth',2.0);
xlabel('$\rho$');ylabel('$w$')
subplot(3,3,6);
plot(vrho,agg_k_vrho,'g','LineWidth',2.0);
xlabel('$\rho$');ylabel('$K$')


subplot(3,3,7);
plot(vsig_y,r_vsig_y,'r','LineWidth',2.0);
xlabel('$\sigma_y$');ylabel('$r$');
subplot(3,3,8);
plot(vsig_y,w_vsig_y,'b','LineWidth',2.0);
xlabel('$\sigma_y$');ylabel('$w$');
subplot(3,3,9);
plot(vsig_y,agg_k_vsig_y,'g','LineWidth',2.0);
xlabel('$\sigma_y$');ylabel('$K$');


saveas(gcf,fig_location+'Fig_Sec3_1_b.png')
end

