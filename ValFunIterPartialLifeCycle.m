function [opt_krt_idx,opt_kwk_idx,opt_cwk,opt_crt] = ...
        ValFunIterPartialLifeCycle(tm_rt,tm_wk,mutilrt,mutilwk,Vterminal,inc,ygrid,survs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%   transitional matrix for work and retire: tm_wk, tm_rt
%   period utility for work adn retire: mutilwk, mutilrt
%   terminal value function (evaluated in the last period no saving): Vterminal
%
% Output:
%   optimal k decision for work and retire: opt_kwk_idx, opt_krt_idx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
global w r beta  ...
    nygrid nwk nrt nkgrid nagrid...
    kgrid ss ...
    iter_max err_max 

iter_max = 1000;
err_max = 10^-6;

iter = 0;
err = 10^6;
Vwk2D = ones(nygrid*nwk,nkgrid);
Vrt2D = ones(nrt-1,nkgrid);

% adjust with survival rate
survs_wk = survs(1:nwk);
msurvs_wk = repmat(survs_wk,[1,nygrid,nkgrid]);
msurvs_wk = reshape(msurvs_wk,[nygrid*nwk,nkgrid]);

survs_rt = survs(nwk+1:nagrid-1);
msurvs_rt = repmat(survs_rt,[1,nkgrid]);

%% Value function iteration
tic
while iter<iter_max && err>err_max
    iter = iter+1;
    
    % rt
    Vrt_temp = tm_rt*[Vrt2D;Vterminal];
    Vrt_temp = Vrt_temp.*msurvs_rt;
    Vrt = mutilrt + beta* permute(repmat(Vrt_temp,[1,1,nkgrid]),[1,3,2]);
    [Vrt_new,opt_krt_idx] = max(Vrt,[],3);
    err1 = max(abs(Vrt_new-Vrt2D),[],'all');
    % wk
    Vwk_temp = tm_wk*[Vwk2D;Vrt_new(1,:)];
    Vwk_temp = Vwk_temp.*msurvs_wk;
    Vwk = mutilwk + beta*permute(repmat(Vwk_temp,[1,1,nkgrid]),[1,3,2]);
    [Vwk_new,opt_kwk_idx] = max(Vwk,[],3);
    err2 = max(abs(Vwk_new-Vwk2D),[],'all');
    
    err = max(err1,err2);
    
    Vwk2D = Vwk_new;
    Vrt2D = Vrt_new;
end

%% Report the performance of VFI
disp('##########################')
if iter<iter_max && err<=err_max
    disp('VFI is done successfully')
else
    if iter>=iter_max
        disp('Exceed rounds of loop')
    end
    if err>err_max
        disp('VF does not converge')
    end

end
toc
disp('##########################')

%% Policy: consumption
opt_kwk_idx_reshape = reshape(opt_kwk_idx,[nygrid,nwk,nkgrid]);
[mywk,minc,mkwk] = ndgrid(ygrid,inc,kgrid);
opt_cwk = w*mywk.*minc + (1+r)*mkwk-kgrid(opt_kwk_idx_reshape);

opt_krt_idx_reshape = reshape(opt_krt_idx,[nrt-1,nkgrid]);
[~,mkrt] = ndgrid(1:nrt-1,kgrid);
opt_crt = (1+r)*mkrt-kgrid(opt_krt_idx_reshape) + ss;

end

