% Parameter definition
r = 0.12;       % Tumor cell growth rate factor
K = 1e6;        % Maximum tumor size 
lambda = 0.693; % Drug half-life of 1 day
Delta = 0.15;   % Drug kill efficiency coefficient

% Initial conditions 
Y0 = [1e4; 5];  % [Initial Tumor Size; Initial Drug Dose]
tspan = [0 30]; % Time span (30 days)

% Analytical Solution for untreated tumor growth using dsolve
fprintf('Running Symbolic dsolve Analysis\n');
syms T_sym(t) r_sym K_sym T0_sym

untreated_ode = diff(T_sym, t) == r_sym * T_sym * (1 - T_sym / K_sym);
cond = T_sym(0) == T0_sym;

analytical_sol = dsolve(untreated_ode, cond);

fprintf('Analytical Solution Formula:\n');
disp(analytical_sol);

% Numerical Solution for the treated system using ode45
[t, Y] = ode45(@(t, Y) tumor_growth_system(t, Y, r, K, Delta, lambda), tspan, Y0);

%Basic Plots
figure;

% Tumor Size Plot
subplot(2, 1, 1);
plot(t, Y(:,1), 'b', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Tumor Size');
title('Tumor Size Over Time');
grid on;

% Drug Concentration Plot
subplot(2, 1, 2);
plot(t, Y(:,2), 'r', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Drug Concentration');
title('Drug Concentration Over Time');
grid on;

% The coupled ODE system function
function dYdt = tumor_growth_system(~, Y, r, K, delta, lambda)
T = Y(1); % The size of the tumor
C = Y(2); % The concentration of the drug used for treatment

dTdt = r*T*(1-(T/K)) - delta*C*T; % ODE calculating tumor growth rate
dCdt = -lambda*C;                 % Explicit negative sign for decay

dYdt = [dTdt; dCdt]; % Combine the derivatives into a single output vector
end