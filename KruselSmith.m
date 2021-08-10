function [coef_new,ols_g_Rsquared,olg_b_Rsquared,business_stat] = ...
    KruselSmith(coef,tm_y,tm_s,p)

global alpha beta sigma delta ...
    infty ygrid kgrid Kgrid sgrid nygrid nkgrid nsgrid nKgrid ...
    iter_max err_max ag ab bg bb rho sig_y...
    Kmin Kmax ...
    obs time_win  pre_burn ygrid K_ss

rng('default');
markov_y = dtmc(tm_y);
y_idx_rand = rand(obs,1);
sim_p0 = round(p*obs);
sim_p0(length(sim_p0)) = obs - sum(sim_p0(1:length(sim_p0)-1));
sim_y_idx = simulate(markov_y,time_win-1,'X0',sim_p0)';
agg_l = sum(sim_y_idx,1);
agg_l0 = agg_l(1);

%% Value Function Iteration
[tm_syK,tm_sK2D] = TransMatrixsyK(tm_s,coef);

[ms,my,mK,mk,mklead] = ndgrid(sgrid,ygrid,Kgrid,kgrid,kgrid);

% Mean K is in logs

mw = (1-alpha)*ms.*(exp(mK)*obs).^alpha*agg_l0^(-alpha);
mr = alpha*ms.*(exp(mK)*obs).^(alpha-1)*agg_l0^(1-alpha)-delta;
mc = mw.*my+(1+mr).*mk-mklead;
mutil = (mc<=0).*(-infty)+(mc>0).*(mc.^(1-sigma)-1)/(1-sigma);
mutil3D = reshape(mutil,[nsgrid*nygrid*nKgrid,nkgrid,nkgrid]);

iter = 0;
err = 10^9;
V2D = ones(nsgrid*nygrid*nKgrid, nkgrid);

while iter<iter_max && err>err_max
    
    iter = iter+1;
    
    Vtemp = tm_syK * V2D;
    V3D = mutil3D + ...
        beta*permute(repmat(Vtemp,[1,1,nkgrid]),[1,3,2]);
    [V2Dnew,opt_k_idx] = max(V3D,[],3);
    err = max(abs(V2Dnew-V2D),[],'all');
    V2D = V2Dnew;
end

opt_k = kgrid(opt_k_idx);
opt_c = nan(size(opt_k));

mc3D = reshape(mc,[nsgrid*nygrid*nKgrid,nkgrid,nkgrid]);

for i = 1:nsgrid*nygrid*nKgrid
    for j = 1:nkgrid
        opt_c(i,j) = mc3D(i,j,opt_k_idx(i,j));
    end
end


%% Simulation

markov_sK = dtmc(tm_sK2D);
sim_k_idx = nan(obs,time_win);
sim_k_idx(:,1) = round(mean([1,nkgrid]));

% define the initial aggregate state
% log of average capital
sim_K0 = log(mean(kgrid(sim_k_idx(:,1))));
sim_K0_idx = find(sim_K0<Kgrid,1,'first');
sim_s0_idx = 1;
sim_sK_idx0 = sub2ind([nsgrid,nKgrid],sim_s0_idx,sim_K0_idx);

% simulate time path of aggregate state

sk_x0 = zeros(nsgrid*nKgrid,1); % 'X0' for starting value of markov simulation
sk_x0(sim_sK_idx0,1)=1;
sim_sK_idx = simulate(markov_sK,time_win-1,'X0',sk_x0)';
sim_sK_idx = repmat(sim_sK_idx,[obs,1]);

[sim_s_idx,sim_K_idx] = ind2sub([nsgrid,nKgrid],sim_sK_idx);
% simulate hetergeneous households

for j = 1:time_win-1
    ind(:,j) = sub2ind([nsgrid,nygrid,nKgrid,nkgrid],...
        sim_s_idx(:,j),sim_y_idx(:,j),sim_K_idx(:,j),sim_k_idx(:,j));
    sim_k_idx(:,j+1) = opt_k_idx(ind(:,j));
end

opt_cvec = reshape(opt_c,[nsgrid*nygrid*nKgrid*nkgrid,1]);
sim_c = opt_cvec(ind);

sim_k = kgrid(sim_k_idx);
sim_s_idx = sim_s_idx(1,:);
sim_K_idx = sim_K_idx(1,:);
sim_K = Kgrid(sim_K_idx);

% change into logs
% mean k is the average capital per person
mean_k = log(sum(sim_k,1)/obs);

mean_k_preburn = mean_k(pre_burn:time_win-1);
mean_klead_preburn = mean_k(pre_burn+1:time_win);
sim_s_idx_preburn = sim_s_idx(pre_burn:time_win-1);

reg_mean_kg = nan(size(sim_s_idx_preburn));
reg_mean_kb = nan(size(sim_s_idx_preburn));
reg_mean_kleadg = nan(size(sim_s_idx_preburn));
reg_mean_kleadb = nan(size(sim_s_idx_preburn));

reg_mean_kg(sim_s_idx_preburn==2) ...
    = mean_k_preburn(sim_s_idx_preburn==2);
reg_mean_kb(sim_s_idx_preburn==1) ...
    = mean_k_preburn(sim_s_idx_preburn==1);
reg_mean_kleadg(sim_s_idx_preburn==2) ...
    = mean_klead_preburn(sim_s_idx_preburn==2);
reg_mean_kleadb(sim_s_idx_preburn==1) ...
    = mean_klead_preburn(sim_s_idx_preburn==1);

ols_g = fitlm(reg_mean_kleadg,reg_mean_kg);
coef_g = table2array(ols_g.Coefficients);
ols_g_Rsquared = ols_g.Rsquared;

ols_b = fitlm(reg_mean_kleadb,reg_mean_kb);
coef_b = table2array(ols_b.Coefficients);
olg_b_Rsquared = ols_b.Rsquared;

% location in the regression table
ag_new = coef_g(1,1);
bg_new = coef_g(2,1);
ab_new = coef_b(1,1);
bb_new = coef_b(2,1);

coef_new = [ag_new, bg_new, ab_new, bb_new];
%% Business Cycle Statistics

% not in logs below
agg_k_pb =  sum(sim_k(:,pre_burn:time_win-1),1);
agg_klead_pb =  sum(sim_k(:,pre_burn+1:time_win),1);
agg_inv_pb = agg_klead_pb - (1-delta)*agg_k_pb;
agg_y_pb = sum(ygrid(sim_y_idx(:,pre_burn+1:time_win)),1);
agg_s_pb = sgrid(sim_s_idx(pre_burn:time_win-1));
agg_output_pb = agg_s_pb.*agg_k_pb.^alpha.*agg_y_pb.^(1-alpha);
agg_c_pb = sum(sim_c,1);
agg_c_pb = agg_c_pb(:,pre_burn:time_win-1);

% change into logs
log_agg_inv = log(agg_inv_pb);
log_agg_y = log(agg_y_pb);
log_agg_s = log(agg_s_pb);
log_agg_output = log(agg_output_pb);
log_agg_c = log(agg_c_pb);

% standard deviation
std_output = std(log_agg_output);
std_inv = std(log_agg_inv);
std_c = std(log_agg_c);
std_s = std(log_agg_s);

% correlation between specific variable and output in logs
corr_output = corrcoef(log_agg_output,log_agg_output);
corr_output = corr_output(1,2);
corr_inv = corrcoef(log_agg_output,log_agg_inv);
corr_inv = corr_inv(1,2);
corr_c = corrcoef(log_agg_output,log_agg_c);
corr_c = corr_c(1,2);
corr_s = corrcoef(log_agg_output,log_agg_s);
corr_s = corr_s(1,2);

business_stat = [std_output,    std_inv,    std_c,      std_s;
        corr_output,    corr_inv,   corr_c,     corr_s];
   
figure();
plot(Kgrid(sim_K_idx));hold on;plot(mean_k);
legend('Perceived Law of Motion','Actual Law of Motion');




end



