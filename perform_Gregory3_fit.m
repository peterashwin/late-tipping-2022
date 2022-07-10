function [lambdas_est, DT_est,ff] = perform_Gregory2_fit(var,par, rolling_window_size)
%% Fit to Gregory fit over rolling windows
%

DT = var.T - par.T0;
R = var.dTdt;

lambdas_est = nan(size(DT));
DT_est = nan(size(DT));
lambdas_se = nan(size(DT));
DT_se = nan(size(DT));
RSq = nan(size(DT));

fit_inds = rolling_window_size:10:length(DT);

for i = fit_inds

    RW = R(i-rolling_window_size+1:i);
    DTW = DT(i-rolling_window_size+1:i);

    % fit to linear model
    mdl=fitlm(DTW, RW);

    tab = mdl.Coefficients;
    coeff = tab.(1);
    coeff_std = tab.(2);
    COV = mdl.CoefficientCovariance;

    lambdaE = coeff(2);
    fE = coeff(1);
    DTE = -fE/lambdaE;

    lambdaE_se = coeff_std(2);
    f_se = coeff_std(1);

    % propogate std uncertainty for q = x/y using formula
    % std_q = |q| sqrt( (std_x/x)^2 + (std_y/y)^2 - 2 * cor_xy/(x*y) )
    DTE_se = abs(DTE) * sqrt( (lambdaE_se/lambdaE).^2 + (f_se/fE).^2 - 2 * COV(1,2) / (fE*lambdaE));

    lambdas_est(i) = lambdaE;
    lambdas_se(i) = lambdaE_se;
    DT_est(i) = DTE;
    DT_se(i) = DTE_se;
    %   R squared
    RSq(i) = mdl.Rsquared.Adjusted;
    
end

%% Plotting
figure();
clf;
f=gcf();
f.Position(3:4)=[330 330];

ts=var.t(fit_inds);

subplot(3,1,1);
plot(ts, lambdas_est(fit_inds), 'k-');
hold on;
plot(ts, lambdas_est(fit_inds)+lambdas_se(fit_inds), 'r-');
plot(ts, lambdas_est(fit_inds)-lambdas_se(fit_inds), 'b-');
plot([0 par.EndTime],[0 0],'k:');
ylim([-0.04, 0.02]);
xlabel('$t$ [year]', 'Interpreter', 'latex');
ylabel('$\lambda$ [$W/m^2/K$]', 'Interpreter', 'latex');
title('Gregory fit', 'Interpreter', 'latex');
hold off;

subplot(3,1,2);
plot(ts, DT_est(fit_inds), 'k-');
hold on;
plot(ts, DT_est(fit_inds)+DT_se(fit_inds), 'r-');
plot(ts, DT_est(fit_inds)-DT_se(fit_inds), 'b-');
plot([0 par.EndTime],[0 0],'k:');
ylim([par.DTminplot par.DTmaxplot]);
xlabel('$t$ [year]', 'Interpreter', 'latex');
ylabel('$\Delta T^*_\mathrm{est}$ [$K$]', 'Interpreter', 'latex');
hold off;

subplot(3,1,3);
plot(var.t(fit_inds), RSq(fit_inds), 'b-');
hold on;
xlabel('$t$ [year]', 'Interpreter', 'latex');
ylabel('$R^2$', 'Interpreter', 'latex');
hold off;

ff=gcf().Number;


end

