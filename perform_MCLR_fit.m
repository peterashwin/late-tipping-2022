function [lambdas, DT_ests, ff] = perform_MCLR_fit(var,par, rolling_window_size)
%perform_MCLR_fit Apply a multicomponent linear regression. That is, fit
%[T,alpha]' = A [T,alpha] + F and derive eigenvalues from the matrix A and
%an estimate for the warming via [T,alpha]_est = - A^{-1} F

%% Obtain warming and albedo change values
DT = var.T - par.T0;
DALB = var.alpha - par.alpha0;
%% Obtain DR = dTdt
DR = var.dTdt;
dALBdt = var.dalphadt;

if sum(isnan(dALBdt)) > 0
    lambdas = nan;
    DT_ests = nan;
   return  
end

%% MC-LR fit over rolling windows

lambdas = nan(length(DT),2);
DT_ests = nan(size(DT));

for i = rolling_window_size:length(DT)
    
   DT_fit = DT(i-rolling_window_size+1:i);
   DALB_fit = DALB(i-rolling_window_size+1:i);
   DR_fit = DR(i-rolling_window_size+1:i);
   dALBdt_fit = dALBdt(i-rolling_window_size+1:i);
   
   X = [ones(size(DT_fit))'; DT_fit'; DALB_fit'];
   Y = [DR_fit'; dALBdt_fit'];
   
   coeff = Y / X;
   
   F = coeff(:,1);
   A = coeff(:,2:3);
   
   X_est = - A \ F;
   lambda = eig(A);
   
   DT_ests(i) = X_est(1);
   lambdas(i,:) = lambda;
    
end

%% Plotting

figure()
f=gcf();
f.Position(3:4)=[330 330];

subplot(2,1,1)
plot(var.t, lambdas, 'r-')
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$\lambda$ [$W/m^2/K$]', 'Interpreter', 'latex')
title('MC-LR')
subplot(2,1,2)
plot(var.t, DT_ests, 'b-')
ylim([-10 100])
xlabel('$t$ [year]', 'Interpreter', 'latex')
ylabel('$\Delta T^*_\mathrm{est}$ [$K$]', 'Interpreter', 'latex')

ff=gcf().Number;



end

