%% script with GEBM run parameters
% March 2022
% define default parameters and then a number of variants as
% pars(1..n)
%

%% Default parameters
% Incoming solar radiation
par.Q0 = 341.3;

% Temperature-dependent equilibrium values for albedo
par.alpha_1 = 0.7;
par.alpha_2 = 0.289;
par.T_alpha = 274.5;
par.K_alpha = 0.1;
par.alpha_0 = @(T,p) p.alpha_1 + (p.alpha_2-p.alpha_1) * ...
    (1 + tanh(p.K_alpha*(T-p.T_alpha)))/2;

% Outgoing radiation
par.sigma = 5.67e-8;

% Temperature-dependent equilibrium values for emissivity
par.eps_1 = 0.5; %0.7; %0.5
par.eps_2 = 0.41; %0.6; %0.41
par.K_eps = 0.5; %0.05; %0.5
par.T_eps = 288;
par.eps_0 = @(T,p) p.eps_1 + (p.eps_2-p.eps_1) * ...
    (1 + tanh(p.K_eps*(T-p.T_eps)))/2;

% Heat capacity 
par.C_T = 5e8;

% seconds per year
par.S = 31556926;

%% Lorenz-63 system
% Standard/original values are sigma=10,beta=8/3 and rho=28
par.sigma_L = 10;
par.rho_L = 28;
par.beta_L = 8/3;
% Initial conditions for Lorenz-63; random ICs in range
par.y0_L = [40 * (rand()-0.5); 50 * (rand()-0.5); 5 + 40 * rand()];

%% Coupling function, strength and timescale for the Lorenz system
par.mu_NV = @(x,p) p.nu_NV * sin(pi*x/20);
par.nu_NV = 2e-2; %5; %5;
par.tau_NV = 6e7; % 3e8 is a timescale of 10 years

%% Time_scale for the albedo evolution
par.tau_alpha = 5e9; %0

%% Equilibrium temperature and albedo
par.Teq = 293;
par.alphaeq = par.alpha_0(par.Teq,par);

%% Initial temp
par.T0 = par.Teq;
par.alpha0=par.alphaeq;

%% Set background CO2 forcing such that Teq is an equilibrium:
par.mu0 = @(p) p.eps_0(p.Teq,p).*p.sigma.*p.Teq.^4-p.Q0 .* (1 - p.alphaeq);

%% Forcing scenario
par.A0 = 5.35;
par.mu = @(t,p) p.A0 * log(4 + 0.* t) + p.mu0(p); % Instantaneous Quadrupling

%% Simulation Setup -- Options for time integration
par.EndTime = 2500;
par.tspan = 0:1:par.EndTime;
par.Name='FastSlow_WARM_4xCO2';

%% Plotting defaults for Delta T
par.DTminplot=-5;
par.DTmaxplot=20;

%%
%
%
% default 
% run WARM 4xCO2
pars(1)=par;

% run COLD 4xCO2
pars(2)=par;
pars(2).Name='FastSlow_COLD_4xCO2';
% change equlibrium Teq
pars(2).Teq = 255;
pars(2).alphaeq = par.alpha_0(pars(2).Teq,pars(2));
% set initial state at eqm
pars(2).T0 = pars(2).Teq;
pars(2).alpha0=pars(2).alphaeq;
pars(2).DTminplot=-10;
pars(2).DTmaxplot=80;

% run COLD 4xCO2 no late tip
pars(3)=pars(2);
pars(3).Name='FastSlow_COLD_4xCO2_slow_Tip';
pars(3).K_eps = 0.1;
pars(3).DTminplot=-10;
pars(3).DTmaxplot=80;

% run COLD 2xCO2
pars(4)=pars(2);
pars(4).Name='FastSlow_COLD_2xCO2';
pars(4).mu = @(t,p) p.A0 * log(2 + 0.* t) + p.mu0(p); % Instantaneous Doubling
pars(4).DTminplot=-5;
pars(4).DTmaxplot=20;




