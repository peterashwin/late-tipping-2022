function [ff] =make_greg_fit_ensemble(var,par, i_first, lambda_bounds, DT_guess)
%make_greg_fit_ensemble: makes a gregory plot (DT,DR) for the raw data, and
%the ensemble-average data, and compute linear regression 'Gregory' fits,
%including error estimates for the estimated equilibrium warming

N = par.EnsembleSize;

%% Obtain temperature and imbalance data

for j=1:N
    DTs(:,j) = var(j).T(i_first:end)-par.T0;
    DRs(:,j) = var(j).dTdt(i_first:end);
end

t = var.t;
t = t(i_first:end);

%%  Perform Gregory fit for each ensemble member individually
%
lambdas = zeros(size(DTs));
lambdas_se = zeros(size(DTs));
DT_ests = zeros(size(DTs));
DT_ests_se = zeros(size(DTs));
RSqs = zeros(size(DTs));

for j = 1:N
    sprintf('Ensemble member= %d',j)
    for i=1:length(t)
        [lambda, lambda_se, DT_est, DT_se, RSq] = Greg_fit(DTs(1:i,j), DRs(1:i,j));
        lambdas(i,j) = lambda;
        lambdas_se(i,j) = lambda_se;
        DT_ests(i,j) = DT_est;
        DT_ests_se(i,j) = DT_se;
        RSqs(i,j) = RSq;
    end
end


% Make figure
figure();
clf;
f=gcf();
f.Position(3:4)=[330 330];

subplot(3,1,1)
plot([t(1),t(end)], [0,0], 'k:', 'linewidth', 2.0,'HandleVisibility','off')
hold on
plot(t, lambdas, '-', 'linewidth', 0.1,'HandleVisibility','off','Color',[0.3 0.3 0.3])
xlabel('$t$ [year]', 'Interpreter', 'latex');
ylabel('$\lambda$ [$W/m^2/K$]', 'Interpreter', 'latex');
ylim(lambda_bounds');
xlim([t(1) t(end)]);

subplot(3,1,2)
plot([t(1),t(end)], [0,0], 'k:', 'linewidth', 2.0,'HandleVisibility','off')
hold on
plot(t,DT_ests, '-', 'linewidth', 0.1,'HandleVisibility','off','Color',[0.3 0.3 0.3])
ylabel('$\Delta T^*_\mathrm{est}$ [$K$]', 'Interpreter', 'latex');
ylim([par.DTminplot par.DTmaxplot]);
xlim([t(1) t(end)]);

subplot(3,1,3)
plot(t, RSqs, '-', 'linewidth', 0.1,'HandleVisibility','off','Color',[0.3 0.3 0.3])
hold on
xlabel('$t$ [year]', 'Interpreter', 'latex');
ylabel('$R^2$', 'Interpreter', 'latex')
ylim([0 1])
xlim([t(1) t(end)])

%% Perform Gregory fit on the ensemble-average
DT_mean = mean(DTs,2);
DR_mean = mean(DRs,2);

lambda_EA = nan(size(t));
lambda_EA_se = nan(size(t));
DT_ests_EA = nan(size(t));
DT_ests_EA_se = nan(size(t));
RSq_EA = nan(size(t));

for i=1:length(t)
    DT_i = DT_mean(1:i);
    DR_i = DR_mean(1:i);
    [lambda, lambda_se, DT_est, DT_se, RSq] = Greg_fit(DT_i, DR_i);
    lambda_EA(i) = lambda;
    lambda_EA_se(i) = lambda_se;
    DT_ests_EA(i) = DT_est;
    DT_ests_EA_se(i) = DT_se;
    RSq_EA(i) = RSq;
end

subplot(3,1,1)
plot(t, lambda_EA, 'r-', 'linewidth', 2.0,'DisplayName','fit to ensemble average')
plot(t, lambda_EA+lambda_EA_se, 'r--','linewidth', 1.0,'HandleVisibility','off')
plot(t, lambda_EA-lambda_EA_se, 'r--', 'linewidth', 1.0,'HandleVisibility','off')

subplot(3,1,2)
plot(t, DT_ests_EA, 'r-', 'linewidth', 2.0,'DisplayName','fit to ensemble average')
plot(t, DT_ests_EA+DT_ests_EA_se, 'r--','linewidth', 1.0,'HandleVisibility','off')
plot(t, DT_ests_EA-DT_ests_EA_se, 'r--', 'linewidth', 1.0,'HandleVisibility','off')

subplot(3,1,3)
plot(t, RSq_EA, 'r-', 'linewidth', 2.0,'DisplayName','fit to ensemble average')

%% Perform Gregory fit on the whole set of data of all ensemble simultaneously

lambda_all = nan(size(t));
lambda_all_se = nan(size(t));
DT_ests_all = nan(size(t));
DT_ests_all_se = nan(size(t));
RSq_all = nan(size(t));

for i=1:length(t)
    DT_all = DTs(1:i,:);
    DR_all = DRs(1:i,:);
    DT_i = DT_all(:);
    DR_i = DR_all(:);
    [lambda, lambda_se, DT_est, DT_se, RSq] = Greg_fit(DT_i, DR_i);
    lambda_all(i) = lambda;
    lambda_all_se(i) = lambda_se;
    DT_ests_all(i) = DT_est;
    DT_ests_all_se(i) = DT_se;
    RSq_all(i) = RSq;
end


%% Plot the theoretical warming

opts1 = optimset('display','off');
F_rhs = @(x) par.Q0.*(1 - par.alpha_0(x,par)) - ...
    par.sigma*par.eps_0(x,par).*x.^4 + par.mu(t(end),par);

DT_eq_real_1 = fsolve( F_rhs , par.T0 + DT_guess(1), opts1) - par.T0;
DT_eq_real_2 = fsolve( F_rhs , par.T0 + DT_guess(2), opts1) - par.T0;

del = 10^(-4);
lambda_eq_real_1 = par.S/par.C_T * ( F_rhs(DT_eq_real_1+par.T0+del)-F_rhs(DT_eq_real_1+par.T0))/norm(del);
lambda_eq_real_2 = par.S/par.C_T * ( F_rhs(DT_eq_real_2+par.T0+del)-F_rhs(DT_eq_real_1+par.T0))/norm(del);

subplot(3,1,1)
plot([t(1), t(end)], lambda_eq_real_1 * [1,1], 'm-', 'linewidth', 2.0,'DisplayName','equilibrium')
plot([t(1), t(end)], lambda_eq_real_2 * [1,1], 'm-', 'linewidth', 2.0,'HandleVisibility','off')

subplot(3,1,2)
plot([t(1),t(end)],DT_eq_real_1 * [1,1], 'm-', 'linewidth', 2.0,'DisplayName','real equilibrium')
plot([t(1),t(end)],DT_eq_real_2 * [1,1], 'm-', 'linewidth', 2.0)

subplot(3,1,3)
plot([t(1), t(end)], -10*[1,1], 'm-', 'linewidth', 2.0, 'DisplayName', 'real equilibrium')

%% Obtain some information of the figure window for saving purposes

cur_fig = gcf();
ff = cur_fig.Number;

subplot(3,1,3)
l=legend;
set(l,'Interpreter','latex','FontSize',5)

end

%% Gregory fitting
function [lambda_est,lambda_se,DT_est,DT_se,RSq] = Greg_fit(DT,DR)
%     X = [ones(size(DT))'; DT'];
%     Y = DR';
%     coeff = Y / X;
%
%     ff = coeff(:,1);
%     lambda = coeff(:,2);
%
%     DT_est = -ff / lambda;

% fit to linear model
mdl=fitlm(DT, DR);

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

lambda_est = lambdaE;
lambda_se = lambdaE_se;
DT_est = DTE;
DT_se = DTE_se;
%   R squared
RSq = mdl.Rsquared.Adjusted;

end
