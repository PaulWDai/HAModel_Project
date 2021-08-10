function err = CalibrationLifeCycle(para,rho,sig_y,inc,con,cons,survs)

beta = para(1);
sigma = para(2);

[mean_c,~,~,~] = ...
    solve_PartialLifeCycle(beta,sigma,rho,sig_y,inc,con,survs);
% starting from 22
mean_c = mean_c(2:length(mean_c));
err = sum((mean_c-cons).^2);

end

