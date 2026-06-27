% ODE solver comparison for the tumor-drug model.
% Compares manual Euler, manual RK4, and built-in ode45 trajectories
% calculates the absolute error tracking over the 30-day window.
% Includes analytical vs numerical comparison, and error and stability
% analysis against closed form solutions for the no-treatment tumour and drug decay cases
clear; clc; close all;

%% parameters
r      = 0.12;    % tumor cell growth rate factor
K      = 1e6;     % maximum tumor size
lambda = 0.693;   % drug half-life of 1 day
Delta  = 0.15;    % drug kill efficiency coefficient
T0     = 1e4;     % initial tumor size
C0     = 5;       % initial drug dose
tspan  = [0 30];  % time span (30 days)
dt     = 0.1;     % fixed step size
t_manual = (tspan(1):dt:tspan(2))';

rhs = @(t, Y) tumor_growth_system(t, Y, r, K, Delta, lambda);
Y0 = [T0; C0];

%% run all three solvers

% ode45
[t_ode, Y_ode] = ode45(rhs, t_manual, Y0);

% Euler Method
[t_euler, Y_euler] = euler_solve(rhs, tspan, dt, Y0);

% RK4
[t_rk4, Y_rk4] = rk4_solve(rhs, tspan, dt, Y0);

%% solver error metrics (vs ode45)
error_euler = abs(Y_ode(:,1) - Y_euler(:,1));
error_rk4   = abs(Y_ode(:,1) - Y_rk4(:,1));

fprintf('\nSolver Error Metrics at Terminal State (Day 30)\n');
fprintf('   Euler Method Global  Error  at Day 30: %.4f cells\n', error_euler(end));
fprintf('   RK4 Method Global Error at Day 30   : %.4f cells\n', error_rk4(end));

%% basic plots
figure;
subplot(2,1,1)
plot(t_ode, Y_ode(:,1), 'LineWidth', 2)
xlabel('Time (days)')
ylabel('Tumour Size (cells)')
title('Tumour Size Over Time')
grid on

subplot(2,1,2)
plot(t_ode, Y_ode(:,2), 'LineWidth', 2)
xlabel('Time (days)')
ylabel('Drug Concentration')
title('Drug Concentration Over Time')
grid on

%% solver comparison — tumour size
figure;
plot(t_ode,   Y_ode(:,1),   'LineWidth', 2); hold on
plot(t_euler, Y_euler(:,1), '--', 'LineWidth', 2)
plot(t_rk4,   Y_rk4(:,1),   ':',  'LineWidth', 2)
xlabel('Time (days)')
ylabel('Tumour Size (cells)')
title('Solver Comparison for Tumour Size')
legend('ode45', 'Euler', 'RK4')
grid on

%% solver comparison - drug concentration
figure;
plot(t_ode,   Y_ode(:,2),   'LineWidth', 2); hold on
plot(t_euler, Y_euler(:,2), '--', 'LineWidth', 2)
plot(t_rk4,   Y_rk4(:,2),   ':',  'LineWidth', 2)
xlabel('Time (days)')
ylabel('Drug Concentration')
title('Solver Comparison for Drug Concentration')
legend('ode45', 'Euler', 'RK4')
grid on

%% error plot (vs ode45)
figure;
plot(t_ode, error_euler, 'LineWidth', 2); hold on
plot(t_ode, error_rk4,   'LineWidth', 2)
xlabel('Time (days)')
ylabel('Absolute Error in Tumour Size (cells)')
title('Tumour Size Error: Euler and RK4 Compared to ode45')
legend('Euler Error', 'RK4 Error')
grid on

%% phase portrait
figure;
plot(Y_ode(:,1), Y_ode(:,2), 'LineWidth', 2)
xlabel('Tumour Size (cells)')
ylabel('Drug Concentration')
title('Phase Portrait: Drug Concentration vs Tumour Size')
grid on

%% parameter sensitivity — initial drug dose
C0_values = [0, 2.5, 5, 10];
figure; hold on;
for i = 1:length(C0_values)
    Y0_test  = [T0; C0_values(i)];
    rhs_test = @(t, Y) tumor_growth_system(t, Y, r, K, Delta, lambda);
    [t_test, Y_test] = ode45(rhs_test, t_manual, Y0_test);
    plot(t_test, Y_test(:,1), 'LineWidth', 2);
end
xlabel('Time (days)'); ylabel('Tumour Size (cells)')
title('Parameter Sensitivity: Effect of Initial Drug Dose C_0')
legend('C_0 = 0','C_0 = 2.5','C_0 = 5','C_0 = 10')
grid on

%% parameter sensivity — drug kill efficiency
Delta_values = [0.05, 0.15, 0.30];
figure; hold on;
for i = 1:length(Delta_values)
    rhs_test = @(t, Y) tumor_growth_system(t, Y, r, K, Delta_values(i), lambda);
    [t_test, Y_test] = ode45(rhs_test, t_manual, [T0; C0]);
    plot(t_test, Y_test(:,1), 'LineWidth', 2);
end
xlabel('Time (days)'); ylabel('Tumour Size (cells)')
title('Parameter Sensitivity: Effect of Drug Kill Efficiency \Delta')
legend('\Delta = 0.05','\Delta = 0.15','\Delta = 0.30')
grid on

%% analytical solutions
C_analysis = C0 * exp(-lambda * t_manual);
T_analysis = (K * T0) ./ (T0 + (K - T0) * exp(-r * t_manual));

% Re-run solvers with no treatment (Delta=0) for fair analytical comparison
rhs_notx = @(t, Y) tumor_growth_system(t, Y, r, K, 0, lambda);
[~, Y_ode_notx]   = ode45(rhs_notx, t_manual, Y0);
[~, Y_euler_notx] = euler_solve(rhs_notx, tspan, dt, Y0);
[~, Y_rk4_notx]   = rk4_solve(rhs_notx, tspan, dt, Y0);

%% analytical vs numerical — tumour (No Treatment) plot
figure;
plot(t_manual, T_analysis,            'k-',  'LineWidth', 2.5); hold on;
plot(t_manual, Y_euler_notx(:,1), 'r--', 'LineWidth', 1.5);
plot(t_manual, Y_rk4_notx(:,1),   'b-.', 'LineWidth', 1.5);
plot(t_manual, Y_ode_notx(:,1),   'g:',  'LineWidth', 2);
legend('Analytical','Euler','RK4','ode45');
xlabel('Time (days)'); 
ylabel('Tumour Size (cells)');
title('Analytical vs Numerical: Tumour Growth (No Treatment)');
grid on;

%% analytical vs numerical — drug decay plot
figure;
plot(t_manual, C_analysis,       'k-',  'LineWidth', 2.5); hold on;
plot(t_manual, Y_euler(:,2), 'r--', 'LineWidth', 1.5);
plot(t_manual, Y_rk4(:,2),   'b-.', 'LineWidth', 1.5);
plot(t_manual, Y_ode(:,2),   'g:',  'LineWidth', 2);
legend('Analytical','Euler','RK4','ode45');
xlabel('Time (days)'); 
ylabel('Drug Concentration');
title('Analytical vs Numerical: Drug Decay');
grid on;

%% error vs time (vs analytical) plot
error_euler_T = abs(Y_euler_notx(:,1) - T_analysis);
error_rk4_T   = abs(Y_rk4_notx(:,1)  - T_analysis);

figure;
plot(t_manual, error_euler_T, 'r-', 'LineWidth', 1.5); hold on;
plot(t_manual, error_rk4_T,   'b-', 'LineWidth', 1.5);
legend('Euler Error','RK4 Error');
xlabel('Time (days)'); 
ylabel('Absolute Error (cells)');
title('Numerical Error vs Analytical Solution');
grid on;

%% step size effect on euler plot
h_vals = [1.0, 0.5, 0.1, 0.05, 0.01];
max_errors = zeros(1, length(h_vals));
for j = 1:length(h_vals)
    hj = h_vals(j);
    tj = (0:hj:30)';
    Tj = zeros(length(tj),1); Tj(1) = T0;
    for i = 1:length(tj)-1
        Tj(i+1) = Tj(i) + hj*(r*Tj(i)*(1 - Tj(i)/K));
    end
    T_anal_j = (K*T0)./(T0+(K-T0)*exp(-r*tj));
    max_errors(j) = max(abs(Tj - T_anal_j));
end
figure;
loglog(h_vals, max_errors, 'wo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Step Size h (days)'); 
ylabel('Max Absolute Error (cells)');
title('Euler Convergence: Max Error vs Step Size');
grid on;

%% stability — large step size plot
h_big = 2.0; t_big = (0:h_big:30)';
T_unstable = zeros(length(t_big),1); T_unstable(1) = T0;
for i = 1:length(t_big)-1
    T_unstable(i+1) = T_unstable(i) + h_big*(r*T_unstable(i)*(1 - T_unstable(i)/K));
end
figure;
plot(t_big,    T_unstable, 'r-o', 'LineWidth', 2); hold on;
plot(t_manual, T_analysis,     'k-',  'LineWidth', 2);
legend('Euler h=2.0','Analytical');
xlabel('Time (days)'); 
ylabel('Tumour Size (cells)');
title('Stability: Euler with Large Step Size');
grid on;

fprintf('Euler  Max Error: %.2f cells | Mean: %.2f cells\n', max(error_euler_T), mean(error_euler_T));
fprintf('RK4    Max Error: %.2f cells | Mean: %.2f cells\n', max(error_rk4_T), mean(error_rk4_T));

%% functions
function dYdt = tumor_growth_system(~, Y, r, K, delta, lambda)
T = Y(1); C = Y(2);
dTdt = r*T*(1-(T/K)) - delta*C*T;
dCdt = -lambda*C;
dYdt = [dTdt; dCdt];
end

function [t, Y] = euler_solve(rhs, tspan, h, Y0)
t = (tspan(1):h:tspan(2))';
Y = zeros(numel(t), numel(Y0));
Y(1,:) = Y0(:)';
for k = 1:numel(t)-1
    dt = t(k+1) - t(k);
    Y(k+1,:) = Y(k,:) + (dt * rhs(t(k), Y(k,:)'))';
end
end

function [t, Y] = rk4_solve(rhs, tspan, h, Y0)
t = (tspan(1):h:tspan(2))';
Y = zeros(numel(t), numel(Y0));
Y(1,:) = Y0(:)';
for k = 1:numel(t)-1
    dt = t(k+1) - t(k);
    yk = Y(k,:)';
    k1 = rhs(t(k),       yk);
    k2 = rhs(t(k)+dt/2,  yk + (dt/2)*k1);
    k3 = rhs(t(k)+dt/2,  yk + (dt/2)*k2);
    k4 = rhs(t(k)+dt,    yk +  dt   *k3);
    Y(k+1,:) = (yk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4))';
end
end
