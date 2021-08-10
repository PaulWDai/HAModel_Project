function [tm_wk,tm_rt] = TransMatrixPartial(tm_wk_y)

global nwk nrt nygrid 

%% Transitional Matrix of Working Age
% y, y', a, a'
tm_wk_y4D = repmat(tm_wk_y,[1,1,nwk,nwk]);
% y, a, y', a'
tm_wk_y4D = permute(tm_wk_y4D,[1,3,2,4]);

tm_wk_a = zeros(nwk,nwk);

for i = 1:nwk-1
    tm_wk_a(i,i+1)=1;
end

% a, a', y, y'
tm_wk_a4D = repmat(tm_wk_a,[1,1,nygrid,nygrid]);
% y, a, y', a'
tm_wk_a4D = permute(tm_wk_a4D,[3,1,4,2]);

tm_wk_4D = tm_wk_y4D.*tm_wk_a4D;

tm_wk = nan(nygrid*nwk,nygrid*nwk+1);
tm_wk(1:nygrid*nwk,1:nygrid*nwk) = reshape(tm_wk_4D,[nygrid*nwk,nygrid*nwk]);

% the last column is the state enter determinsitically into retiring age
% threshold
tm_wk(1:nygrid*nwk,nygrid*nwk+1) = 1-sum(tm_wk(1:nygrid*nwk,1:nygrid*nwk),2);

tm_wk = sparse(tm_wk);

%% Transitional Matrix of Retiring Age
tm_rt = zeros(nrt-1,nrt);
for i = 1:nrt-1
    tm_rt(i,i+1)=1;
end


tm_wk = sparse(tm_wk);
tm_rt = sparse(tm_rt);
end

