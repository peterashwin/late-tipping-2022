function ff=make_time_series_figures_ensemble(var,par)
%make_time_series_figures makes time series figures for the variables

N=par.EnsembleSize;
%% Temperature
for j=1:N
    DTs(:,j) = var(j).T - par.T0;
end
t=var.t;

%% create figure
figure();
clf;
f=gcf();
f.Position(3:4)=[330 330];

%plot(t,DTs(:,ceil( rand(1,10)*length(DTs(1,:)) )), 'k-', 'linewidth', 0.1)
plot(t,DTs(:,:), 'k-', 'linewidth', 0.1)
hold on
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$\Delta T$ [K]', 'Interpreter', 'latex')
plot(t, mean(DTs,2), 'r-', 'linewidth', 2.0)


%% Obtain the expected (average) equilibrium temperature
opts1 = optimset('display','off');
DT_eq_real = fsolve( @(x) par.Q0.*(1 - par.alpha_0(x,par)) - ...
    par.sigma*par.eps_0(x,par).*x.^4 + par.mu(t(end),par), par.T0, opts1) - par.T0;

%plot([t(1) t(end)], DT_eq_real*[1;1], 'b--')

cur_fig = gcf();
ff = cur_fig.Number;

end

