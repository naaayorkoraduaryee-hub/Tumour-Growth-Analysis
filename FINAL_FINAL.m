% Applied Programming Project — Final Optimized Code
clear; clc;

%% parameter definition
r = 0.12;       % Tumor cell growth rate factor
K = 1e6;        % Maximum tumor size
lambda = 0.693; % Drug half-life of 1 day
Delta = 0.15;   % Drug kill efficiency coefficient

%% simulation settings
T0 = 1e4;             % Initial tumor size
C0 = 5;                % Drug dose per administration
tspan = [0 30];        % Time span (30 days)
dt = 0.1;              % Time step size for RK4

params.r = r;
params.K = K;
params.lambda = lambda;
params.Delta = Delta;

%% solver comparison (Euler vs RK4 vs ode45)
t_manual = (tspan(1):dt:tspan(2))';
rhs = @(t, Y) tumor_growth_system(t, Y, params);
Y0 = [T0; C0];

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

%% control and toxicity definitions
% controlled = tumor size at day 30 does not exceed its starting size T0
% toxicity = area under the drug concentration curve
is_controlled = @(Y_final) Y_final(1) <= T0;

%% optimization: find the maximum (least frequent) dosing interval that still controls the tumor
optimal_interval = find_optimal_interval(params, T0, C0, tspan, dt, 1, 15, 0.05);

[t_opt, Y_opt] = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, optimal_interval, params);
toxicity_opt = trapz(t_opt, Y_opt(:,2));

[t_safe, Y_safe] = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, 2, params);
toxicity_safe = trapz(t_safe, Y_safe(:,2));

fprintf('\nOptimization Results\n');
fprintf('Minimum effective (maximum tolerable) dose interval: %.2f days\n', optimal_interval);
fprintf('Final tumor size at optimal interval: %.2f cells (T0 = %.0f)\n', Y_opt(end,1), T0);
fprintf('Drug exposure (AUC) at optimal interval: %.2f\n', toxicity_opt);
fprintf('Drug exposure (AUC) at a 2-day comparison interval: %.2f\n', toxicity_safe);
fprintf('Toxicity reduction vs 2-day dosing: %.1f%%\n', 100*(1 - toxicity_opt/toxicity_safe));

%% Level 1A: solution vs time and state-variable evolution
figure('Name','Optimal Treatment State Evolution');
try theme(gcf,'light'); catch, end

subplot(2,1,1);
plot(t_opt, Y_opt(:,1), 'LineWidth', 1.8);
hold on;
yline(T0, '--k', 'Initial tumour size');
xlabel('Time (days)');
ylabel('Tumour cells, T(t)');
title(sprintf('Tumour evolution at optimal interval = %.2f days', optimal_interval));
grid on;

subplot(2,1,2);
plot(t_opt, Y_opt(:,2), 'LineWidth', 1.8);
xlabel('Time (days)');
ylabel('Drug concentration, C(t)');
title('Drug concentration evolution with repeated dosing');
grid on;

set(gcf,'Color','w');
saveas(gcf,'fig2_solution_state_evolution.png');

%% Level 1B: analytical vs numerical comparison
% For the drug equation dC/dt = -lambda*C, the exact single-dose solution is
% C(t) = C0*exp(-lambda*t). Compare it with Euler, RK4 and ode45.
C_exact = C0 .* exp(-params.lambda .* t_manual);

figure('Name','Analytical vs Numerical Drug Concentration');
try theme(gcf,'light'); catch, end
plot(t_manual, C_exact, 'k-', 'LineWidth', 2.2);
hold on;
plot(t_euler, Y_euler(:,2), '--', 'LineWidth', 1.4);
plot(t_rk4, Y_rk4(:,2), '-.', 'LineWidth', 1.4);
plot(t_ode, Y_ode(:,2), ':', 'LineWidth', 1.8);
xlabel('Time (days)');
ylabel('Drug concentration, C(t)');
title('Analytical vs numerical solution for drug decay');
legend('Analytical','Euler','RK4','ode45','Location','northeast');
grid on;
set(gcf,'Color','w');
saveas(gcf,'fig3_analytical_vs_numerical.png');

%% Level 2A: solver comparison plots
figure('Name','Solver Comparison');
try theme(gcf,'light'); catch, end

subplot(2,1,1);
plot(t_ode, Y_ode(:,1), 'k-', 'LineWidth', 2.0);
hold on;
plot(t_euler, Y_euler(:,1), '--', 'LineWidth', 1.3);
plot(t_rk4, Y_rk4(:,1), '-.', 'LineWidth', 1.3);
xlabel('Time (days)');
ylabel('Tumour cells, T(t)');
title('Tumour solution: Euler vs RK4 vs ode45');
legend('ode45 reference','Euler','RK4','Location','best');
grid on;

subplot(2,1,2);
plot(t_ode, Y_ode(:,2), 'k-', 'LineWidth', 2.0);
hold on;
plot(t_euler, Y_euler(:,2), '--', 'LineWidth', 1.3);
plot(t_rk4, Y_rk4(:,2), '-.', 'LineWidth', 1.3);
xlabel('Time (days)');
ylabel('Drug concentration, C(t)');
title('Drug solution: Euler vs RK4 vs ode45');
legend('ode45 reference','Euler','RK4','Location','best');
grid on;

set(gcf,'Color','w');
saveas(gcf,'fig4_solver_comparison.png');

%% Level 2B: error plots
error_euler_C = abs(Y_ode(:,2) - Y_euler(:,2));
error_rk4_C   = abs(Y_ode(:,2) - Y_rk4(:,2));

figure('Name','Solver Absolute Errors');
try theme(gcf,'light'); catch, end

subplot(2,1,1);
semilogy(t_manual, max(error_euler, eps), '--', 'LineWidth', 1.5);
hold on;
semilogy(t_manual, max(error_rk4, eps), '-.', 'LineWidth', 1.5);
xlabel('Time (days)');
ylabel('|Tumour error|');
title('Absolute tumour error relative to ode45');
legend('Euler error','RK4 error','Location','best');
grid on;

subplot(2,1,2);
semilogy(t_manual, max(error_euler_C, eps), '--', 'LineWidth', 1.5);
hold on;
semilogy(t_manual, max(error_rk4_C, eps), '-.', 'LineWidth', 1.5);
xlabel('Time (days)');
ylabel('|Drug error|');
title('Absolute drug-concentration error relative to ode45');
legend('Euler error','RK4 error','Location','best');
grid on;

set(gcf,'Color','w');
saveas(gcf,'fig5_error_plots.png');

%% Level 2C: phase portrait
figure('Name','Tumour Drug Phase Portrait');
try theme(gcf,'light'); catch, end
plot(Y_opt(:,2), Y_opt(:,1), 'LineWidth', 1.8);
hold on;
plot(Y_opt(1,2), Y_opt(1,1), 'go', 'MarkerFaceColor','g', 'MarkerSize',7);
plot(Y_opt(end,2), Y_opt(end,1), 'ro', 'MarkerFaceColor','r', 'MarkerSize',7);
xlabel('Drug concentration, C(t)');
ylabel('Tumour cells, T(t)');
title(sprintf('Phase portrait at %.2f-day dosing interval', optimal_interval));
legend('Trajectory','Start','End','Location','best');
grid on;
set(gcf,'Color','w');
saveas(gcf,'fig6_phase_portrait.png');

%% Level 3: animated dynamic time-evolving simulation
% Set CREATE_ANIMATION to false when you only want the static figures.
CREATE_ANIMATION = true;
if CREATE_ANIMATION
    animation_file = 'fig7_dynamic_treatment_simulation.gif';
    fig_anim = figure('Name','Dynamic Treatment Simulation');
    try theme(fig_anim,'light'); catch, end
    set(fig_anim,'Color','w');

    frame_step = max(1, round(0.25/dt));  % approximately one frame per 0.25 day
    frame_indices = unique([1:frame_step:length(t_opt), length(t_opt)]);

    for frame_number = 1:length(frame_indices)
        k = frame_indices(frame_number);
        clf(fig_anim);

        subplot(2,1,1);
        plot(t_opt(1:k), Y_opt(1:k,1), 'LineWidth', 1.8);
        hold on;
        plot(t_opt(k), Y_opt(k,1), 'o', 'MarkerFaceColor','auto');
        yline(T0, '--k', 'Initial tumour size');
        xlim(tspan);
        ylim([0, 1.05*max(Y_opt(:,1))]);
        xlabel('Time (days)');
        ylabel('Tumour cells, T(t)');
        title(sprintf('Dynamic treatment response: day %.1f', t_opt(k)));
        grid on;

        subplot(2,1,2);
        plot(t_opt(1:k), Y_opt(1:k,2), 'LineWidth', 1.8);
        hold on;
        plot(t_opt(k), Y_opt(k,2), 'o', 'MarkerFaceColor','auto');
        xlim(tspan);
        ylim([0, 1.05*max(Y_opt(:,2))]);
        xlabel('Time (days)');
        ylabel('Drug concentration, C(t)');
        title(sprintf('Repeated dose interval = %.2f days', optimal_interval));
        grid on;

        drawnow;
        frame = getframe(fig_anim);
        [image_data, colour_map] = rgb2ind(frame2im(frame), 256);
        if frame_number == 1
            imwrite(image_data, colour_map, animation_file, 'gif', ...
                'LoopCount', Inf, 'DelayTime', 0.08);
        else
            imwrite(image_data, colour_map, animation_file, 'gif', ...
                'WriteMode', 'append', 'DelayTime', 0.08);
        end
    end
end

%% Simulation results and parameter Sensitivity analysis
%% Section 1: dose-response sweep
intervals = 1:0.5:15;
final_tumor = zeros(length(intervals),1);
peak_tumor  = zeros(length(intervals),1);
auc_drug    = zeros(length(intervals),1);
controlled  = zeros(length(intervals),1);

for i = 1:length(intervals)
    [t_i, Y_i] = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, intervals(i), params);
    final_tumor(i) = Y_i(end,1);
    peak_tumor(i)  = max(Y_i(:,1));
    auc_drug(i)    = trapz(t_i, Y_i(:,2));
    controlled(i)  = is_controlled(Y_i(end,:));
end

disp('Interval | FinalTumor | PeakTumor | DrugAUC | Controlled');
disp([intervals', final_tumor, peak_tumor, auc_drug, controlled]);

figure("Name","Dose Response");
try theme(gcf,'light'); catch, end   
subplot(2,1,1);
plot(intervals, final_tumor, '-o', 'LineWidth', 1.5);
hold on;
yline(T0, '--k', 'Starting tumour size');
xline(optimal_interval, ':k', sprintf('Optimal = %.2f days', optimal_interval));
xlabel('Dosing interval (days)');
ylabel('Final tumour size');
title('Final tumour size vs. dosing interval'); grid on;

subplot(2,1,2);
plot(intervals, auc_drug, '-s', 'LineWidth', 1.5, 'Color', [0.8 0.3 0.1]);
xlabel('Dosing interval (days)');
ylabel('Drug exposure');
title('Drug exposure vs. dosing interval'); grid on;
set(gcf,'Color','w');
saveas(gcf, 'fig1_dose_response.png');

%% Section 2: parameter sensitivity
pt_changes = [-20, -10, 10, 20];
sens_r = zeros(length(pt_changes),4); sens_K = zeros(length(pt_changes),4);
sens_lambda = zeros(length(pt_changes),4); sens_Delta = zeros(length(pt_changes),4);

for i = 1:length(pt_changes)
    pt = pt_changes(i);

    p = params; p.r = params.r * (1 + pt/100);
    [oi, ft, tox] = run_one_scenario(p, T0, C0, tspan, dt);
    sens_r(i,:) = [pt, oi, ft, tox];

    p = params; p.K = params.K * (1 + pt/100);
    [oi, ft, tox] = run_one_scenario(p, T0, C0, tspan, dt);
    sens_K(i,:) = [pt, oi, ft, tox];

    p = params; p.lambda = params.lambda * (1 + pt/100);
    [oi, ft, tox] = run_one_scenario(p, T0, C0, tspan, dt);
    sens_lambda(i,:) = [pt, oi, ft, tox];

    p = params; p.Delta = params.Delta * (1 + pt/100);
    [oi, ft, tox] = run_one_scenario(p, T0, C0, tspan, dt);
    sens_Delta(i,:) = [pt, oi, ft, tox];
end

param_names = {'r','K','lambda','Delta'};
sens_all    = {sens_r, sens_K, sens_lambda, sens_Delta};

row10 = find(pt_changes == 10);
fprintf('\nElasticity at +10%% (OptInterval, Toxicity):\n');
for p_idx = 1:4
    data = sens_all{p_idx};
    el_oi  = (100*(data(row10,2)-optimal_interval)/optimal_interval) / data(row10,1);
    el_tox = (100*(data(row10,4)-toxicity_opt)/toxicity_opt) / data(row10,1);
    fprintf('%-7s %6.3f  %6.3f\n', param_names{p_idx}, el_oi, el_tox);
end

figure('Name','Sensitivity of Optimal Interval');
try theme(gcf,'light'); catch, end
hold on;
for p_idx = 1:4
    data = sens_all{p_idx};
    plot([data(1,2), data(4,2)], [p_idx, p_idx], '-o', 'LineWidth', 2, 'MarkerFaceColor', 'w');
end
xline(optimal_interval, '--k', sprintf('Baseline = %.2f days', optimal_interval));
set(gca, 'YTick', 1:4, 'YTickLabel', param_names);
ylim([0.5, 4.5]);
xlabel('Optimal dosing interval (days)');
title('Sensitivity of optimal interval'); grid on;
set(gcf,'Color','w');


%% Section 2b: Continuous optimal-interval sensitivity sweep

delta_range = linspace(0.05, 0.30, 26);   
opt_interval_vs_delta = zeros(size(delta_range));
for i = 1:length(delta_range)
    p = params;
    p.Delta = delta_range(i);
    opt_interval_vs_delta(i) = find_optimal_interval(p, T0, C0, tspan, dt, 1, 15, 0.05);
end

lambda_range = linspace(0.20, 1.40, 26);  
opt_interval_vs_lambda = zeros(size(lambda_range));
for i = 1:length(lambda_range)
    p = params;
    p.lambda = lambda_range(i);
    opt_interval_vs_lambda(i) = find_optimal_interval(p, T0, C0, tspan, dt, 1, 15, 0.05);
end


figure('Name', 'Continuous Sensitivity Sweep');
try theme(gcf,'light'); catch, end

subplot(1,2,1);
plot(delta_range, opt_interval_vs_delta, '-o');
hold on;
plot(params.Delta, optimal_interval, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Kill efficiency, \Delta');
ylabel('Optimal dosing interval (days)');
title('Optimal Interval vs \Delta');
legend('Swept result', 'Baseline case', 'Location', 'southeast');
grid on;

subplot(1,2,2);
plot(lambda_range, opt_interval_vs_lambda, '-o', 'Color', [0.8 0.3 0.1]);
hold on;
plot(params.lambda, optimal_interval, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Drug decay rate, \lambda');
ylabel('Optimal dosing interval (days)');
title('Optimal Interval vs \lambda');
legend('Swept result', 'Baseline case', 'Location', 'southwest');
grid on;

set(gcf,'Color','w');


%% Section 3: r vs Delta interaction grid
r_values = linspace(params.r * 0.6, params.r * 1.4, 7);
delta_values = linspace(params.Delta * 0.6, params.Delta * 1.4, 7);
heat = zeros(length(delta_values), length(r_values));

for a = 1:length(delta_values)
    for b = 1:length(r_values)
        p = params; p.r = r_values(b); p.Delta = delta_values(a);
        [~, Y_g] = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, optimal_interval, p);
        heat(a,b) = Y_g(end,1);
    end
end

figure('Name','r vs Delta Grid');
try theme(gcf,'light'); catch, end
imagesc(r_values, delta_values, heat);
set(gca, 'YDir', 'normal'); colorbar; hold on;
plot(params.r, params.Delta, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Growth rate, r');
ylabel('Kill efficiency, Delta');
title(sprintf('Final tumour size, interval = %.2f days', optimal_interval));
legend('Baseline (r, Delta)');
set(gcf,'Color','w');


%% FUNCTIONS
function dYdt = tumor_growth_system(~, Y, params)
    T = Y(1);
    C = Y(2);
    dTdt = params.r * T * (1 - (T / params.K)) - params.Delta * C * T;
    dCdt = -params.lambda * C;
    dYdt = [dTdt; dCdt];
end

function [t, Y] = rk4_multidose_solve(rhs, tspan, h, T0, dose_amount, interval, params)
    t = (tspan(1):h:tspan(2))';
    N = length(t);
    Y = zeros(N, 2);
    Y(1, :) = [T0, dose_amount];
    dosing_times = tspan(1):interval:tspan(2);

    for k = 1:(N - 1)
        t_curr = t(k);
        yk = Y(k, :)';

        if k > 1 && any(abs(t_curr - dosing_times) < (h / 2))
            yk(2) = yk(2) + dose_amount;
        end

        k1 = rhs(t_curr,       yk,               params);
        k2 = rhs(t_curr + h/2, yk + (h/2) * k1, params);
        k3 = rhs(t_curr + h/2, yk + (h/2) * k2, params);
        k4 = rhs(t_curr + h,   yk + h * k3,     params);

        Y(k+1, :) = (yk + (h/6) * (k1 + 2*k2 + 2*k3 + k4))';
    end
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

function optimal_interval = find_optimal_interval(params, T0, C0, tspan, dt, low, high, tol)
    is_controlled = @(Y_final) Y_final(1) <= T0;
    [~, Y_low]  = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, low, params);
    [~, Y_high] = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, high, params);
    if ~is_controlled(Y_low(end,:))
        optimal_interval = NaN; return;
    end
    if is_controlled(Y_high(end,:))
        optimal_interval = high; return;
    end
    while (high - low) > tol
        mid = (low + high) / 2;
        [~, Y_mid] = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, mid, params);
        if is_controlled(Y_mid(end,:))
            low = mid;
        else
            high = mid;
        end
    end
    optimal_interval = low;
end

function [opt_interval, final_tumor, toxicity] = run_one_scenario(params, T0, C0, tspan, dt)
    opt_interval = find_optimal_interval(params, T0, C0, tspan, dt, 1, 15, 0.05);
    if isnan(opt_interval)
        final_tumor = NaN; toxicity = NaN; return;
    end
    [t, Y] = rk4_multidose_solve(@tumor_growth_system, tspan, dt, T0, C0, opt_interval, params);
    final_tumor = Y(end,1);
    toxicity = trapz(t, Y(:,2));
end