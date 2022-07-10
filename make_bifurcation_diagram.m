function fn=make_bifurcation_diagram(var,par)
%make_bifurcation_diagram Makes a bifurcation diagram for the parameter
%values used in the simulation. Because of the equation this is just a plot
%of the function mu = eps(T) sigma T^4 - Q0 ( 1 - alpha_0(T) )
%Also plots the evolution of the simulation projected onto this diagram.
%The stars indicate the starting and end points

%% Obtain the graph (mus, Ts) via the function for mu
Ts = linspace(200, 350);

mus = par.eps_0(Ts,par) .* par.sigma .* Ts.^4 - par.Q0 * (1 - par.alpha_0(Ts,par));

%% Make figure
figure()
plot(mus,Ts)
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('$T_*$', 'Interpreter', 'latex')

%% Add current time series to the plot
T = var.T;
mu = par.mu(var.t,par);

hold on
plot(mu, T, 'r:')
plot(par.mu0(par), par.T0, 'b*')
plot(mu(end), T(end), 'r*')

fn=gcf().Number;

end

