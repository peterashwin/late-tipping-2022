%% make_all_figures_ensemble
% This script analyses the output of the numerical simulations made with
% the script doGEBMsimulationbatch.m.

%% Close figures
close all

%% make filename to read
file_name = [path name];
load([file_name '.mat'])

%% Make time series figures
f1=make_time_series_figures_ensemble(var,par);
savefigure(name,f1,printfigs);

%% Make histograms
f2=make_histogram_heat_map_figures_ensemble(var,par);
savefigure(name,f2,printfigs);

%% Make Gregory plot & Gregory fit

f3=make_greg_fit_ensemble(var,par, i_first, lambda_bounds, DT_guess);
savefigure(name,f3,printfigs);
