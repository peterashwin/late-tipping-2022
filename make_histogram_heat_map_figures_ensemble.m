function ff=make_histogram_heat_map_figures_ensemble(var,par)
%make_histogram_heat_map_figures_ensemble: makes histograms for all times
%of the spread in temperature warming putting them in the bins as dictated
%by the array edges.

N=par.EnsembleSize;

edges_DT=par.DTminplot:par.DTstep:par.DTmaxplot;

%% Temperature
for j=1:N
    Ts(:,j)=var(j).T;
end
t=var.t;

%% Cut the first element (i.e. at t0)
% this is to prevent a large spike at the initial condition
del_elems = 10;

t = t(del_elems:end);
Ts = Ts(del_elems:end,:);

%% Compute the warming
DTs = Ts - par.T0;

%% Histograms at various times

edges{2} = edges_DT;
edges{1} = t;

%% create figure
figure();
clf;
f=gcf();
f.Position(3:4)=[330 330];

% We need to convert data to a single array for hist3 to work
DT_values = DTs(:);
t_times = repmat(t,length(Ts(1,:)),1);
X = [t_times, DT_values];
hist3(X, 'CdataMode', 'auto', 'edges', edges );
ylabel('$\Delta T$ [K]', 'Interpreter', 'latex');
xlabel('$t$ [year]', 'Interpreter', 'latex');
shading interp;
axis tight;
a = colorbar;
a.Label.String = 'Number of occurences';
a.Label.Interpreter = 'latex';
view(2);


%% Obtain the expected (average) equilibrium temperature
opts1 = optimset('display','off');
DT_eq_real = fsolve( @(x) par.Q0.*(1 - par.alpha_0(x,par)) - ...
    par.sigma*par.eps_0(x,par).*x.^4 + par.mu(t(end),par), par.T0, opts1) - par.T0;

cur_fig = gcf();
ff = cur_fig.Number;

end

