function [mutilrt,mutilwk,Vterminal] = PeriodUtilityPartialLifeCycle(ygrid,inc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Calculation of period utility and the termination value (no asset
%   left) at given state. Notice: the life cycle is defined as work age and
%   retirement (no labor productivity shock, reduce dimensionality).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global w r sigma ...
    nygrid nwk nrt nkgrid...
     artgrid kgrid ss infty


%   wk = work;  rt = retired;
%   a = age, k = capital this period, klead = capital next

%% Work
[mywk,minc,mkwk,mkleadwk] = ndgrid(ygrid,inc,kgrid,kgrid);
mcwk = w*mywk.*minc + (1+r)*mkwk-mkleadwk;
mutilwk = (mcwk<=0)*(-infty) + (mcwk>0).*mcwk.^(1-sigma)/(1-sigma);
mutilwk = reshape(mutilwk,[nygrid*nwk,nkgrid,nkgrid]);
%% Retire
[~,mkrt,mkleadrt] = ndgrid(artgrid(1:nrt-1),kgrid,kgrid);
mcrt = (1+r)*mkrt-mkleadrt + ss;
mutilrt = (mcrt<=0)*(-infty) + (mcrt>0).*mcrt.^(1-sigma)/(1-sigma);
%   no asset holding to next period at end of the life
cterminal = (1+r)*kgrid+ss;
Vterminal = (cterminal<=0)*(-infty) + (cterminal>0).*cterminal.^(1-sigma)/(1-sigma);
end

