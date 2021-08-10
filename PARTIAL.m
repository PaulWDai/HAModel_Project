% Project 2
% Author: Weifeng Dai, ID: 1700018613

clear;
clc;
close;
%% Parameters

% Benchmark:
beta = 1/(1+0.04);
sigma = 3;
rho = .9;
sig_y = .4;

sigma_vec = [1.001,5];
rho_new = 0;
sig_y_vec = [.3,.5];

amin = 21;          % grid for age
abar = 65;          % the last age of work
amax = 80;

nagrid = amax-amin+1;
nwk = abar-amin+1;      % number of age of work
nrt = amax-(abar+1)+1;  % number of age of retired
agrid = linspace(amin,amax,nagrid);

% read files of life cycle data;
incfile = fopen('incprofile.txt','r');
consfile = fopen('consprofile.txt','r');
survsfile = fopen('survs.txt');
formatSpec = '%f';
inc = fscanf(incfile,formatSpec);
con = fscanf(consfile,formatSpec);
con = reshape(con,[2,.5*length(con)])';

for i = 22:87
    cons(1,i-21) = mean(con(4*(i-22)+1:4*(i-22)+4,2));
end
cons = cons/cons(1,1);
cons = cons(1,1:nagrid-1);

survs = fscanf(survsfile,formatSpec);
survs = survs(length(survs)-nagrid+1:length(survs),1);

% settings
fig_location = "/Users/weifengdai/Dropbox/My Mac (MacBook Pro)/Documents/Undergrad/2021Spring/Heterogeneous Macro/Project/Project2/figure/";
%% Benchmark
[mean_c_benchmark,mean_k_benchmark,~,~] = ...
        solve_PartialLifeCycle(beta,sigma,rho,sig_y,inc,con,survs);

    
%% Calibration

%{
para0 = [beta,sigma];
cal_para = fminsearch(@(para)CalibrationLifeCycle(para,rho,sig_y,inc,con,cons,survs),para0);

[mean_c,~,~,~] = ...
    solve_PartialLifeCycle(beta,sigma,rho,sig_y,inc,con,survs);
% starting from 22
mean_c = mean_c(2:length(mean_c));
err = sum(abs(mean_c-cons));

%}

%% Variation of Sigma
for idx_sigma = 1:length(sigma_vec)
    [mean_c,mean_k,~,~] = ...
        solve_PartialLifeCycle(beta,sigma_vec(idx_sigma),rho,sig_y,inc,con,survs);
    mean_c_sigma(idx_sigma,:) = mean_c;
    mean_k_sigma(idx_sigma,:) = mean_k;
    %std_c_sigma(idx_sigma,:) = std_c;
    %std_k_sigma(idx_sigma,:) = std_k;
end

%% Variation of Rho

[mean_c,mean_k,~,~] = ...
    solve_PartialLifeCycle(beta,sigma,rho_new,sig_y,inc,con,survs);
mean_c_rho = mean_c;
mean_k_rho = mean_k;
%std_c_rho = std_c;
%std_k_rho = std_k;


%% Variation of Sig_y

for idx_sig_y = 1:length(sig_y_vec)
    [mean_c,mean_k,~,~] = ...
        solve_PartialLifeCycle(beta,sigma,rho,sig_y_vec(idx_sig_y),inc,con,survs);
    mean_c_sig_y(idx_sig_y,:) = mean_c;
    mean_k_sig_y(idx_sig_y,:) = mean_k;
    %std_c_sig_y(idx_sig_y,:) = std_c;
    %std_k_sig_y(idx_sig_y,:) = std_k;
end

%% Variation of Survival Rate

survs_new = ones(nagrid,1);
survs_new(length(survs_new),1)=0;
[mean_c,mean_k,std_c,std_k] = ...
    solve_PartialLifeCycle(beta,sigma,rho,sig_y,inc,con,survs_new);
mean_c_survs = mean_c;
mean_k_survs = mean_k;
std_c_survs = std_c;
std_k_survs = std_k;

%% Life Cycle Consumption and Asset
figure();
sgtitle('Life Cycle Consumption')
subplot(2,2,1);
plot(agrid,mean_c_sigma(1,:),'r:','LineWidth',3.0);hold on;
plot(agrid,mean_c_benchmark,'b-','LineWidth',3.0);hold on;
plot(agrid,mean_c_sigma(2,:),'g-.','LineWidth',3.0);hold on;
title('Variation of $\sigma$');
legend('$\sigma=1$','$\sigma=3$','$\sigma=5$','location','Northwest');
xlabel('Age');ylim([0.8 3]);

subplot(2,2,2);
plot(agrid,mean_c_rho,'r:','LineWidth',3.0);hold on;
plot(agrid,mean_c_benchmark,'b-','LineWidth',3.0);hold on;
title('Variation of $\rho$');
legend('$\rho=0$','$\rho=0.9$','location','Northwest');
xlabel('Age');ylim([0.8 3]);

subplot(2,2,3);
plot(agrid,mean_c_sig_y(1,:),'r:','LineWidth',3.0);hold on;
plot(agrid,mean_c_benchmark,'b-','LineWidth',3.0);hold on;
plot(agrid,mean_c_sig_y(2,:),'g-.','LineWidth',3.0);hold on;
title('Variation of $\sigma_y$');
legend('$\sigma_y=0.3$','$\sigma_y=0.4$','$\sigma_y = 0.5$','location','Northwest');
xlabel('Age');ylim([0.8 3]);

subplot(2,2,4);
plot(agrid,mean_c_benchmark,'r:','LineWidth',3.0);hold on;
plot(agrid,mean_c_survs,'b-','LineWidth',3.0);hold on;
title('Variation of Survival Rate');
legend('With Survival Rate','Without Survival Rate','location','Northwest');
xlabel('Age');ylim([0.8 3]);

saveas(gcf,fig_location+'Fig_Sec2_Consumption.png');

figure();
sgtitle('Life Cycle Asset')
subplot(2,2,1);
plot(agrid,mean_k_sigma(1,:),'r:','LineWidth',3.0);hold on;
plot(agrid,mean_k_benchmark,'b-','LineWidth',3.0);hold on;
plot(agrid,mean_k_sigma(2,:),'g-.','LineWidth',3.0);hold on;
title('Variation of $\sigma$');
legend('$\sigma=1$','$\sigma=3$','$\sigma=5$','location','Northwest');
xlabel('Age');ylim([0 25]);

subplot(2,2,2);
plot(agrid,mean_k_rho,'r:','LineWidth',3.0);hold on;
plot(agrid,mean_k_benchmark,'b-','LineWidth',3.0);hold on;
title('Variation of $\rho$');
legend('$\rho=0$','$\rho=0.9$','location','Northwest');
xlabel('Age');ylim([0 25]);

subplot(2,2,3);
plot(agrid,mean_k_sig_y(1,:),'r:','LineWidth',3.0);hold on;
plot(agrid,mean_k_benchmark,'b-','LineWidth',3.0);hold on;
plot(agrid,mean_k_sig_y(2,:),'g-.','LineWidth',3.0);hold on;
title('Variation of $\sigma_y$');
legend('$\sigma_y=0.3$','$\sigma_y=0.4$','$\sigma_y = 0.5$','location','Northwest');
xlabel('Age');ylim([0 25]);

subplot(2,2,4);
plot(agrid,mean_k_benchmark,'r:','LineWidth',3.0);hold on;
plot(agrid,mean_k_survs,'b-','LineWidth',3.0);hold on;
title('Variation of Survival Rate');
legend('With Survival Rate','Without Survival Rate','location','Northwest');
xlabel('Age');ylim([0 25]);

saveas(gcf,fig_location+'Fig_Sec2_Asset.png');

