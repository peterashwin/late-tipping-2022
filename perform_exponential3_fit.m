function [lambdas_est, DT_est,ff] = perform_exponential_fit(var,par, rolling_window_size)
% Fit to exponential decay

%% Obtain warming and albedo change values
DT = var.T - par.T0;

%% Fit to decaying exponential fit over rolling windows

fitopts=statset('MaxIter',300,'Display','iter','TolFun',1e-12);

lambdas_est = nan(size(DT));
DT_est = nan(size(DT));
lambdas_se = nan(size(DT));
DT_se = nan(size(DT));
Rsq = nan(size(DT));

bestim=nan(3,length(DT));
bse=nan(3,length(DT));

scale_b3 = 1/rolling_window_size;
fit_inds = rolling_window_size:10:length(DT);

for i = fit_inds

    X = var.t(i-rolling_window_size+1:i);
    X = (X - X(1))/rolling_window_size;

    Y = DT(i-rolling_window_size+1:i);

    % Model to fit to
    modelfun = @(b,x) (b(1)+b(2)*exp(-b(3)*x(:,1)));

    % Initial guess
    beta0 = [ Y(end) -1.0 1.0];

    mdl = fitnlm(X, Y, modelfun, beta0, 'Options', fitopts);

    tab = mdl.Coefficients;
    bestim(:,i)=tab.(1);
    bse(:,i)=tab.(2);

    RS=mdl.Rsquared.Adjusted;

    lambdas_est(i) = -1*bestim(3,i)/rolling_window_size;
    lambdas_se(i) = bse(3,i)/rolling_window_size;
    DT_est(i) = bestim(1,i);
    DT_se(i) = bse(1,i);
    RSq(i) = RS;
    
end

%bestim(3,:) = bestim(3,:)*scale_b3;
%bse(3,:) = bse(3,:)*scale_b3;


%% Plotting

figure()
f=gcf();
f.Position(3:4)=[330 330];

ts=var.t(fit_inds);

subplot(3,1,1)
hold on
plot(ts, lambdas_est(fit_inds), 'k-');
plot(ts, lambdas_est(fit_inds)+lambdas_se(fit_inds), 'r-');
plot(ts, lambdas_est(fit_inds)-lambdas_se(fit_inds), 'b-');
plot([0 par.EndTime],[0 0],'k:');
ylim([-0.04, 0.02])
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$\lambda$ [$W/m^2/K$]', 'Interpreter', 'latex')
title('Fit to decaying exponential', 'Interpreter', 'latex')
hold off

subplot(3,1,2)
hold on
plot(ts, DT_est(fit_inds), 'k-');
plot(ts, DT_est(fit_inds)+DT_se(fit_inds), 'r-');
plot(ts, DT_est(fit_inds)-DT_se(fit_inds), 'b-');
plot([0 par.EndTime],[0 0],'k:');
ylim([par.DTminplot par.DTmaxplot]);
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$\Delta T^*_\mathrm{est}$ [$K$]', 'Interpreter', 'latex')
hold off

subplot(3,1,3);
plot(var.t(fit_inds), RSq(fit_inds), 'b-');
hold on;
xlabel('$t$ [year]', 'Interpreter', 'latex');
ylabel('$R^2$', 'Interpreter', 'latex');
hold off;

ff=gcf().Number;


end

