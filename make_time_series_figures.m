function [ff1,ff2]=make_time_series_figures(var, par)
%make_time_series_figures makes time series figures for the variables

%% Temperature and albedo
DT = var.T - par.T0;
Dalpha = var.alpha - par.alpha0;

figure();
f=gcf();
f.Position(3:4)=[700 330];
subplot(2, 2, 1);
plot(var.t, DT, 'r-','LineWidth',1.5)
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$\Delta T$ [K]', 'Interpreter', 'latex')

subplot(2, 2, 3);
plot(var.t, Dalpha, 'b-','LineWidth',1.5)
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$\Delta \alpha$', 'Interpreter', 'latex')

subplot(2, 2,[2 4]);
plot(DT, Dalpha, 'c-','LineWidth',1.5)
title('Temperature and Albedo', 'Interpreter', 'latex')
xlabel('$\Delta T$ [K]', 'Interpreter', 'latex')
ylabel('$\Delta \alpha$', 'Interpreter', 'latex')

% also plot the graphs alpha = alpha_0(T) and alpha = (1+mu/Q0) -
% eps(T)sigmaT^4/Q0
Ts = floor(min(var.T)-10):0.1:ceil(max(var.T)+20);
alpha_alphas = par.alpha_0(Ts,par);
alpha_Ts = (1 + par.mu(var.t(end),par)/par.Q0) - ...
    par.eps_0(Ts,par).*par.sigma.*Ts.^4/par.Q0;
hold on
plot(Ts-par.T0, alpha_alphas-par.alpha0, 'r:','LineWidth',1.5)
plot(Ts-par.T0, alpha_Ts-par.alpha0, 'b:','LineWidth',1.5)

xlim([par.DTminplot par.DTmaxplot]);
%xlim([-10 80]);

legend('orbit', 'nullcline alpha', 'nullcline temperature', 'Location','southwest')

ff1=gcf().Number;

%% Lorenz System
figure()
subplot(3,4,1)
plot(var.t, var.x_L, 'r-')
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$x$', 'Interpreter', 'latex')

subplot(3,4,5)
plot(var.t, var.y_L, 'r-')
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')

subplot(3,4,9)
plot(var.t, var.z_L, 'r-')
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$z$', 'Interpreter', 'latex')

subplot(3,4,[2 3 4 6 7 8 10 11 12])
scatter3(var.x_L, var.y_L, var.z_L)
xlabel('$xL$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$z$', 'Interpreter', 'latex')
title('Lorenz variables', 'Interpreter', 'latex')

ff2=gcf().Number;

end

