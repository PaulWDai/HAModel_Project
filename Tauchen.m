function [tm_final,ygrid] = Tauchen(nstate,rho,sig_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: 
%   Discrete the AR1 process using Tauchen's method

% Input:
%   nstate: number of state
%   rho:    autocorr.
%   sig_y:  std. of transitory component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logygrid = linspace(-3*sig_y,3*sig_y,nstate);
logystep = (max(logygrid)-min(logygrid))/(nstate-1);
tm = nan(nstate,nstate);
for state1 = 1:nstate
    for state2 = 1:nstate
        ystate2 = logygrid(state2);
        ystate1 = logygrid(state1);
        temp = ystate2-rho*ystate1;
        if state2>1 && state2<nstate
            tm(state1,state2) = ...
                normcdf((temp+logystep/2)/sig_y) ...
                - normcdf((temp-logystep/2)/sig_y);
        else
            if state2==1
                tm(state1,state2) = ...
                    normcdf((temp+logystep/2)/sig_y);
            
            else
                tm(state1,state2) = ...
                    1-normcdf((temp-logystep/2)/sig_y);
            end
        end
    end
end

tm_final = tm;

ygrid = exp(logygrid);

end

