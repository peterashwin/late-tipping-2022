function [t,y] = heun_integration(ode, tspan, y0, options)
%heun_integration Performs time integration of ode via heun integration
% ode must be a function of y and t @(t,y) ode(t,y)
% ode must be a function handle
% tspan is the vector of desired output times
% y0 is the initial condition


%% Noise and time step information
etas = options.eta; % amplitude of noise
dt = options.dt; % timestep

n = size(y0,1);

% Make the noise in the right dimensions depending on amount of variables
% all other dimensions will be set to 0
if ( options.is_fastslow )
    eta = etas(1:2);
else
    eta = etas(1);
end
if ( options.is_chaotic )
    eta = [eta; etas(end-2:end)];
end


%% Actual intergration

t0 = tspan(1);
k = 2; % Index for which value in tspan we are looking for

Tout = t0;
Yout = y0';

Yp = y0;
tp = tspan(1);


tn = t0 + dt;
while tn < tspan(end)
    
    % Compute noise
    Wn = sqrt(dt) * randn(n,1);
    
    %RK2/heun steps
    m1=feval(ode,tp,Yp);
    Yt=Yp+m1*dt+eta.*Wn;
    m2=feval(ode,tn,Yt);
    Yn=Yp+0.5*(m1+m2)*dt+eta.*Wn;

    
    % Add data to the to-be-outputted data array
    % but try to find values for the actual values in tspan based on linear
    % interpolation of the current and previous values
    while ( k <= length(tspan) ) && ( tn >= tspan(k) )
        tk = tspan(k);
        yk = Yp + (Yn-Yp)/(tn-tp) * (tk - tp);
        
        Tout = [Tout;tn];
        Yout = [Yout;yk'];
        
        k = k + 1; 
    end

    %update previous data
    tp=tn;
    Yp=Yn;
    
    % Update for the time
    tn = tn + dt;
end

%% Make the output
t = Tout;
y = Yout;

end
