function p_final = StationaryStatePDF(tm_y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Find the stationary distibution of state of a markov chain

% Output:
%   p_final: the stationary distibution of state of a markov chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global iter_max err_max nygrid
iter_markov = 0;
err_markov = 10^9;

p = ones(1,nygrid)*1/nygrid;

while iter_markov < iter_max && err_markov>err_max
    iter_markov = iter_markov+1;
    pnew = p* tm_y;
    err_markov = max(abs(p-pnew));
    p = pnew;
end

p_final = p;

end

