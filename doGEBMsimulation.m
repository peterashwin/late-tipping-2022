%% doGEBMsimulation
%
% This script runs a numerical simulation for the following energy balance
% model, which potentially has fast-slow behaviour and/or coupling to a
% chaotic Lorenz system.
% C_T dT/dt = Q0 (1 - alpha) - epsilon(T) sigma T^4 + mu + mu_NV
% tau_alpha dalpha/dt = alpha_0(T) - alpha
% mu_NV = nu_NV *sin(pi x/20), where x is the first component of a
% Lorenz-63 model, i.e.
% tau_NV dx/dt = sigma_L (y-x)
% tau_NV dy/dt = x (rho_L - z) - y
% tau_NV dz/dt = xy - beta_L z
%

%% Start with a clean slate
close all
clear all
%% Start timer
tic
%% Model Parameters -- global energy balance model load parameters

run('GEBMrunparams.m');

for i=1:length(pars)
    %% run through each of the parameter values
    par=pars(i);

    %% Initial state vector
    par.y0 = [par.T0; par.alpha0];

    %% Simulation Setup -- which numerical integration method to use
    options.time_integrator = 'ode45';

    % The following options are possible
    % ode45: use MATLAB's ode45 function
    % Options for ode45
    options.ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); % for ode45

    % heun: heun's method for solving odes (can include noise)
    % Options for heun
    % options.time_integrator = 'heun';
    % NOTE: the time step dt should be small enough! This is increasingly
    % important as the value tau_NV is smaller as then very small time steps
    % are needed to properly simulate the Lorenz system.

    options.eta=[1e-6;1e-6;0;0;0]; % amplitude of noise (columns: T, alpha, x_L, y_L, z_L)
    options.dt=0.1; % time step (in years!) of the heun integrator

    %% Call the function that runs the actual numerical simulation
    [var] = GEBMsimulator(par,options);

    %% Saving
    path = '../Data/';
    name = par.Name;
    file_name = [path name];

    save( [file_name '.mat'], 'par', 'var', 'options');

    %% Stop timer
    toc

end