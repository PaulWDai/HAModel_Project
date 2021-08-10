function [opt_k_idx_final,opt_c_final,V2D_final] = ...
    ValueFunInfinite(r,w,rho,sigma,sig_y,tau,ind_T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: 
%   Solving the value function and the policy function of infinitely-lived
%   agents optimizations porblem

% Input:
%   r,w: rental price of captial and wage
%   rho, sigma, sig_y: parameters of autocorr. of income, CRRA coef., std.
%   of transitory income shock
%   tau: capital income tax
%   ind_T: individual level of transfer from gov.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global alpha beta delta ...
    infty ygrid kgrid nygrid nkgrid ...
    iter_max err_max tm_y

[my,mk,mklead] = ndgrid(ygrid, kgrid, kgrid);

mc = w* my +(1+(1-tau)*r)*mk-mklead + ind_T;
mutil = (mc>0).*mc.^(1-sigma)/(1-sigma) + ...
    (mc<=0).*(-infty);

iter = 0;
err = 10^9;

V2D = ones(nygrid, nkgrid);

while iter <iter_max && err>err_max
    iter = iter+1;
    Vtemp = tm_y * V2D;
    V3D = mutil + ...
        beta*permute(repmat(Vtemp,[1,1,nkgrid]),[1,3,2]);
    [V2Dnew,opt_k_idx] = max(V3D,[],3);
    err = max(abs(V2Dnew-V2D),[],'all');
    V2D = V2Dnew;
end

V2D_final = V2D;
opt_k_idx_final = opt_k_idx;


for i = 1:nygrid
    for j = 1:nkgrid
        opt_c(i,j) = mc(i,j,opt_k_idx(i,j));
    end
end

opt_c_final = opt_c;
end

