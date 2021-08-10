function []= ...
    Fig_Sec3_1_a(k_density_markov,k_density_monte_carlo,sample_k,sample_c)

global kgrid fig_location

figure();
subplot(1,3,1);title('Comparison: Asset Distribution')
plot(kgrid,k_density_markov,'b.','LineWidth',1.0);hold on;
plot(kgrid,k_density_monte_carlo,'r-','LineWidth',1.0);
legend('Markov','Monte Carlo');...
    xlabel('$k$','Interpreter','latex'), ylabel('PDF');title('wealth distribution');ylim([0 0.03]);

subplot(1,3,2);
histogram(sample_k);title('Asset Inequality');

subplot(1,3,3);
histogram(sample_c);title('Consumption Inequality');

saveas(gcf,fig_location+'Fig_Sec3_1_a.png');
end

