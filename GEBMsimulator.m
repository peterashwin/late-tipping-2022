function [var] = GEBMsimulator(par,options)
%GEBMsimulator. Performs numerical time integration of energy balance
%model.

    %% Construct the right initial condition
    % if the chaotic Lorenz system is needed then we need that in our
    % initial conditions, otherwise we do not need that.
    % Similarly, only consider alpha as seprate variable if tau_alpha is
    % not 0.
    
    is_chaotic = (par.nu_NV ~= 0);
    is_fastslow = (par.tau_alpha ~= 0);
    options.is_chaotic = is_chaotic;
    options.is_fastslow = is_fastslow;
    
    if ( is_fastslow )
        y0 = par.y0(1:2);
    else
        y0 = par.y0(1);
    end
    
    if ( is_chaotic )
        % just add the initial conditions for the lorenz system to it
        y0 = [y0; par.y0_L];
    end
        

    %% Selection of the time integrator plus time integration
    
    switch options.time_integrator
        case 'ode45'
            [t,y] = ode45(@(t,y) balance_model(t,y,par), par.tspan, y0, options.ode_opts);
        case 'heun'
            [t,y] = heun_integration(@(t,y) balance_model(t,y,par), par.tspan, y0, options);     
        otherwise
            warning("Unexpected numerical time integration method chosen that has not been programmed in.")
    end


   %% Obtain the desired outputs
   
   % Obtain the (instantaneous) time derivatives via the ODE function
   dydt = nan(size(y));
   for i=1:length(t)
       dydt(i,:) = balance_model( t(i), y(i,:), par);
   end
   % Then obtain the separate variables
  
   T = y(:,1);
   dTdt = dydt(:,1);
   % if tau_alpha is 0 then alpha is not a state variable
   if ( is_fastslow )
      alpha = y(:,2);
      dalphadt = dydt(:,2);
   else
      alpha = par.alpha_0(T,par);
      dalphadt = nan(size(T));
   end
   % If nu_NV = 0 then we do not have the Lorenz system part
   if (is_chaotic)
       x_L = y(:,end-2);
       y_L = y(:,end-1);
       z_L = y(:,end);
   else
       x_L = nan(size(T));
       y_L = nan(size(T));
       z_L = nan(size(T));
   end
   
   
    %% Put everything together
    var.t = t;
    var.T = T;
    var.dTdt = dTdt;
    var.alpha = alpha;
    var.dalphadt = dalphadt;
    var.x_L = x_L;
    var.y_L = y_L;
    var.z_L = z_L;

end


%% Create the right-hand side of the equation
% Here there are a few scenarios
% If tau_alpha = 0 then we need to set alpha = alpha_0(T) manually and need
% to eliminate alpha as dynamic system component
% If nu_NV then we do not need the Lorenz-63 model as it does not influence
% the temperature evolution
function dydt = balance_model(t, y, par)
    
    is_chaotic = (par.nu_NV ~= 0);
    is_fastslow = (par.tau_alpha ~= 0);
    
    % Temperature evolution and albedo evolution    
    dTdt = @(T,alpha,mu) par.Q0 * (1-alpha) - ...
        par.eps_0(T,par) .* par.sigma .* T.^4 + mu;
    
    if is_fastslow
        T = y(1);
        alpha = y(2);
        
        dydt = par.S * [ dTdt(T,alpha, par.mu(t,par)) / par.C_T; ...
            (par.alpha_0(T,par)-alpha) / par.tau_alpha ];
    else
        T = y(1);
        
        dydt = par.S * dTdt(T,par.alpha_0(T,par),par.mu(t,par)) / par.C_T;
    end
    
    
    % Lorenz-system if needed
    if is_chaotic
       z = y(end-2:end); % obtain lorenz variables
       % The evolution of the Lorenz system:
       dzdt = [ par.sigma_L * ( z(2) - z(1) ); ...
           z(1) .* ( par.rho_L - z(3) ) - z(2) ; ...
           z(1) .* z(2) - par.beta_L * z(3) ] * par.S / par.tau_NV;
       % Coupling to the energy balance model
       mu_NV = par.mu_NV(z(1),par);
       % add this to dydt's first component
       dydt(1) = dydt(1) + mu_NV * par.S/par.C_T;
       % Combine everything into one vector for dydt
       dydt = [ dydt; dzdt];
    end
    
end
    