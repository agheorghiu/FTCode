f = fit(dp(1:end)', results2(1:end)', 'gauss2', 'Robust', 'on');

%errs = 0.01 * ones(1, 12);%errs(1) = 0;

xpoints = 0.05 : 0.005 : 0.15;
plot(dp, results1, xpoints, f(xpoints), 'linewidth', 1.3);
hold on;
plot(dp, results2, 'ko', 'markers', 5);
errorbar(dp,results2, errs, 'ko', 'markers', 6);

xlabel('Probability of depolarization error');
ylabel('Probability of acceptance');
legend('Unencoded state', 'Encoded state', 'Data points Steane code', 'Location', 'northeast');
%set(findall(gca, 'Type', 'Line'),'LineWidth',1.3);
set(gca, 'XTick', 0.05 : 0.01 : 0.14);
xlim([0.05 0.14]);