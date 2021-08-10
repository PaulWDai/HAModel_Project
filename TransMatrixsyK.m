function [tm_syK,tm_sK2D] = TransMatrixsyK(tm_s,coef)
global  ...
     Kgrid  nygrid  nsgrid nKgrid...
     tm_y ag ab bg bb Kmin Kmax 


ag = coef(1);
bg = coef(2);
ab = coef(3);
bb = coef(4);

Kleadg = ag + bg*Kgrid;
Kleadb = ab + bb*Kgrid;

Kleadgidx = nan(1,nKgrid);
Kleadbidx = nan(1,nKgrid);

for i = 1:nKgrid
    idx_g_temp = find(Kleadg(1,i)<Kgrid,1);
    % Consider the index of Klead exceeds the max of Kgrid
    if isempty(idx_g_temp)==1
        Kleadgidx(1,i) = nKgrid;
    else
        Kleadgidx(1,i) = idx_g_temp;
    end
    
    idx_b_temp = find(Kleadb(1,i)<Kgrid,1);
    if isempty(idx_b_temp)==1
        Kleadbidx(1,i) = nKgrid;
    else
        Kleadbidx(1,i) = idx_b_temp;
    end
end

tm_Kg = zeros(nKgrid,nKgrid);
tm_Kb = zeros(nKgrid,nKgrid);

for iter_K = 1:nKgrid
    tm_Kg(iter_K,Kleadgidx(1,iter_K)) = 1;
    tm_Kb(iter_K,Kleadbidx(1,iter_K)) = 1;
end

% s, s', K, K'
tm_sK = nan(nsgrid,nsgrid,nKgrid,nKgrid);
tm_sK(1,1,:,:) = tm_s(1,1)*tm_Kb; % bad state to bad state
tm_sK(1,2,:,:) = tm_s(1,2)*tm_Kb; % bad state to good state
tm_sK(2,1,:,:) = tm_s(2,1)*tm_Kg; % good state to bad state
tm_sK(2,2,:,:) = tm_s(2,2)*tm_Kg; % good state to good state

% s, K, s', K'
tm_sK2D = permute(tm_sK,[1,3,2,4]);
tm_sK2D = reshape(tm_sK2D,[nsgrid*nKgrid,nsgrid*nKgrid]);

% s, s', K, K', y, y'
tm_sK6D = repmat(tm_sK,[1,1,1,1,nygrid,nygrid]);
% s, y, K, s', y', K'
tm_sK6D = permute(tm_sK6D,[1,5,3,2,6,4]);

% y, y', s, s', K, K'
tm_y6D = repmat(tm_y,[1,1,nsgrid,nsgrid,nKgrid,nKgrid]);
% s, y, K, s', y', K'
tm_y6D = permute(tm_y6D,[3,1,5,4,2,6]);

tm6D = tm_sK6D.*tm_y6D;
% (s,y,K), (s',y', K')
tm_syK = reshape(tm6D,[nsgrid*nygrid*nKgrid,nsgrid*nygrid*nKgrid]);

end

