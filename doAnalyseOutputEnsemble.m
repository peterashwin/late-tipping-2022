%% doAnalyseOutputEnsemble
% This script analyses the output of the numerical simulations made with
% the script doGEBMsimulationbatch.m.
% and prints figures using make_all_figures_ensemble.m

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
path = '../Data/';

%% uncomment lines that you want to run
%
% name = 'EnsembleAbrupt4xCO2_T293';
% i_first = 1;
% DT_guess = [0;0];
% lambda_bounds = [-0.2,0.1];
% run('make_all_figures_ensemble.m');

%%
name = 'EnsembleAbrupt4xCO2_always';
i_first = 76;
DT_guess = [70;70];
lambda_bounds = [-0.2,0.1];
run('make_all_figures_ensemble.m');

%%
% name = 'EnsembleAbrupt4xCO2_to75y';
% i_first = 76;
% DT_guess = [0;70];
% lambda_bounds = [-0.2,0.1];
% run('make_all_figures_ensemble.m');
