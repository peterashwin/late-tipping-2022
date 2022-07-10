%% doAnalyseOutput
% This script analyses the output of the numerical simulations made with
% the script doGEBMsimulation.m.
% and prints figures using make_all_figures.m

%% Start with a clean slate
clear all

%% OPTIONAL: TURN OFF WARNINGS
warning('off','all')
warning
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%%
printfigs=true;

%% Obtain data from file
path = '..\Data\';

%% Uncomment lines that you want to run
%
name = 'FastSlow_WARM_4xCO2';
run('make_all_figures.m');

name = 'FastSlow_COLD_2xCO2';
run('make_all_figures.m');

name = 'FastSlow_COLD_4xCO2';
run('make_all_figures.m');

name = 'FastSlow_COLD_4xCO2_slow_Tip';
run('make_all_figures.m');
