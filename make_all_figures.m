%% make_all_figures.m
% This script analyses the output of the numerical simulations made with
% the script doGEBMsimulation.m.
% It uses that data to do the following
% 1. Make Bifurcation Diagram with mu as bifurcation parameter
% 2. Make time series figures and plot orbits
% 3. Apply linear 'Gregory' fit to data to estimate equilibrium warming
% 4. Apply MC-LR fit to data (See [Bastiaansen et al, 2020]) to estimate
% equilibrium warming
% 5. Apply a nonlinear fit of a decaying exponential to estimate
% equilibrium warming

% close figures
close all

%% make filename to read
file_name = [path name];
load([file_name '.mat'])

%% Construct bifurcation diagram
f1=make_bifurcation_diagram(var,par);
%savefigure(name,f1,printfigs);

%% Make time series figures
[f2,f3]=make_time_series_figures(var, par);
savefigure(name,f2,printfigs);
%savefigure(name,f3,printfigs);

%% Gregory Fits
rolling_fit_window = 150;
[lambdas_GREG, DT_ests_GREG,f4] = perform_Gregory3_fit(var,par, rolling_fit_window);
savefigure(name,f4,printfigs);

% %% MC-LR Fits
% rolling_fit_window = 150;
% [lambdas_MCLR, DT_ests_MCLR,f5] = perform_MCLR_fit(var,par, rolling_fit_window);
% savefigure(name,f5,printfigs);

%% Exponential fit
rolling_fit_window = 150;
[lambdas_expfit, DT_ests_expfit,f6] = perform_exponential3_fit(var,par, rolling_fit_window);
savefigure(name,f6,printfigs);
