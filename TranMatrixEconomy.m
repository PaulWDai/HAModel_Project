function sim_tm2D = TranMatrixEconomy(opt_k_idx,tm_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: 
%   Find oout the transitional matrix of entire economy, given the
%   exogenous transitional matrix of state, i.e., tm_y and the policy
%   function opt_k_idx

% Input:
%   opt_k_idx: policy function (row is state other than k, col is k state,
%   value is the next period k index).
%   tm_y: transitional matrix exogenously given. e.g. income process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nygrid nkgrid

sim_tm = nan(nygrid,nkgrid,nygrid,nkgrid);

for idx_y = 1:nygrid
    for idx_k = 1:nkgrid
        for idx_ylead = 1:nygrid
            for idx_klead = 1:nkgrid
                sim_tm(idx_y,idx_k,idx_ylead,idx_klead) = ...
                    tm_y(idx_y,idx_ylead)*(opt_k_idx(idx_y,idx_k)==idx_klead);
            end
        end
    end
end

sim_tm2D = reshape(sim_tm,[nygrid*nkgrid,nygrid*nkgrid]);

end

