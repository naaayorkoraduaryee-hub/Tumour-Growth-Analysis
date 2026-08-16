clear; clc; close all;

%% Graph 1 Chinese Hamster V79 fibroblast tumor 
t_data = [3.46 4.58 5.67 6.64 7.63 8.41 9.32 10.27 11.19 12.39 13.42 15.19 16.24 17.23 18.18 19.29 21.23 21.99 24.33 25.58 26.43 27.44 28.43 30.49 31.34 32.34 33.00 35.20 36.34 37.29 38.50 39.67 41.37 42.58 45.39 46.38 48.29 49.24 50.19 51.14 52.10 54.00 56.33 57.33 59.38];

V_data = [0.0158 0.0264 0.0326 0.0445 0.0646 0.0933 0.1454 0.2183 0.2842 0.4977 0.6033 0.8441 1.2163 1.4470 2.3298 2.5342 3.0064 3.4044 3.2046 4.5241 4.3459 5.1374 5.5376 4.8946 5.0660 6.1494 6.8548 5.9668 6.6945 6.6395 6.8971 7.2966 7.2268 6.8815 8.0993 7.2112 7.0694 7.4971 6.9974 6.7219 7.0523 7.1095 7.0694 8.0562 7.2268];

V0 = V_data(1);
t0 = t_data(1);

%Fit data 
cost = @(p) nmse_cost(p, t_data, V_data, t0, V0);

params0 = [0.3, 7];   % initial guess
opts = optimset('TolX', 1e-10, 'TolFun', 1e-12, 'MaxFunEvals', 20000, 'MaxIter', 20000);
params_opt = fminsearch(cost, params0, opts);

a = params_opt(1);
b = params_opt(2);
NMSE = cost(params_opt);

fprintf('Fitted parameters:\n  a = %.5f\n  b = %.4f\n  NMSE = %.6f\n', a, b, NMSE);

%Solve ODE
tspan = linspace(t0, 60, 600);
f = @(t, V) a .* V .* (1 - V ./ b);
[t_sol, V_sol] = ode45(f, tspan, V0);

%Plot 
figure();
plot(t_data, V_data, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold on;
plot(t_sol, V_sol, 'r-', 'LineWidth', 1.5);
hold off;

xlabel('Days');
ylabel('Tumor Volume (mm^3)');
title('Logistic Chinese Hamster V79 fibroblast tumor ');
xlim([0 60]);
ylim([0 8]);
box on;

annotation('textbox', [0.3 0.01 0.4 0.08], 'String',sprintf('NMSE = %.6f\na = %.5f\nb = %.4f', NMSE, a, b),'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);



%% Graph 2 (Logistic Glioblastoma)
t_data2 = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70];
V_data2 = [148 175 225 326 426 559 593 680 792 849 927 999 1085 1138 1217];

V0_2 = V_data2(1);
t0_2 = t_data2(1);

cost2 = @(p) nmse_cost(p, t_data2, V_data2, t0_2, V0_2);

params0_2 = [0.06, 1300];
opts2 = optimset('TolX', 1e-10, 'TolFun', 1e-12, 'MaxFunEvals', 20000, 'MaxIter', 20000);
params_opt2 = fminsearch(cost2, params0_2, opts2);

a2 = params_opt2(1);
b2 = params_opt2(2);
NMSE2 = cost2(params_opt2);

fprintf('\nSecond dataset fitted parameters:\n  a = %.6f\n  b = %.4f\n  NMSE = %.7f\n', a2, b2, NMSE2);

tspan2 = linspace(t0_2, 70, 600);
f2 = @(t, V) a2 .* V .* (1 - V ./ b2);
[t_sol2, V_sol2] = ode45(f2, tspan2, V0_2);

figure();
plot(t_data2, V_data2, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold on;
plot(t_sol2, V_sol2, 'r-', 'LineWidth', 1.5);
hold off;

xlabel('Days');
ylabel('Tumor Volume (mm^3)');
title('Logistic Glioblastoma');
xlim([0 70]);
ylim([0 1300]);
box on;

annotation('textbox', [0.3 0.01 0.4 0.08], 'String',sprintf('NMSE = %.7f\na = %.6f\nb = %.4f', NMSE2, a2, b2),'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);


%% Graph 3 (NOD/SCID Mice)
t_data3 = [0.00 8.95 22.93 26.01 33.11 36.89 39.86 43.92 47.97];
V_data3 = [2.52 26.52 122.65 143.96 197.00 297.02 333.38 563.22 819.57];

V0_3 = V_data3(1);
t0_3 = t_data3(1);

cost3 = @(p) nmse_cost(p, t_data3, V_data3, t0_3, V0_3);

params0_3 = [0.145, 989];
opts3 = optimset('TolX', 1e-10, 'TolFun', 1e-12, 'MaxFunEvals', 20000, 'MaxIter', 20000);
params_opt3 = fminsearch(cost3, params0_3, opts3);

a3 = params_opt3(1);
b3 = params_opt3(2);
NMSE3 = cost3(params_opt3);

fprintf('\nThird dataset fitted parameters:\n  a = %.6f\n  b = %.4f\n  NMSE = %.7f\n', a3, b3, NMSE3);

tspan3 = linspace(t0_3, 48, 600);
f3 = @(t, V) a3 .* V .* (1 - V ./ b3);
[t_sol3, V_sol3] = ode45(f3, tspan3, V0_3);

figure();
plot(t_data3, V_data3, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold on;
plot(t_sol3, V_sol3, 'r-', 'LineWidth', 1.5);
hold off;

xlabel('Days');
ylabel('Tumor Volume (mm^3)');
title('Logistic NOD/SCID Mice');
xlim([0 50]);
ylim([0 900]);
box on;

annotation('textbox', [0.3 0.01 0.4 0.08], 'String',sprintf('NMSE = %.7f\na = %.5f\nb = %.2f', NMSE3, a3, b3),'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);


%% Graph 4 (BALB/Slc-nu/nu Mice)

t_data4 = [6.93 9.08 10.97 12.87 14.92 16.91 19.98 24.90 27.97 30.90 35.97];
V_data4 = [177.96 219.14 291.20 450.02 594.13 682.37 826.48 1073.54 1217.66 1747.06 2082.35];

V0_4 = V_data4(1);
t0_4 = t_data4(1);

cost4 = @(p) nmse_cost(p, t_data4, V_data4, t0_4, V0_4);

params0_4 = [0.136, 2494];
opts4 = optimset('TolX', 1e-10, 'TolFun', 1e-12, 'MaxFunEvals', 20000, 'MaxIter', 20000);
params_opt4 = fminsearch(cost4, params0_4, opts4);

a4 = params_opt4(1);
b4 = params_opt4(2);
NMSE4 = cost4(params_opt4);

fprintf('\nFourth dataset fitted parameters:\n  a = %.6f\n  b = %.4f\n  NMSE = %.7f\n', a4, b4, NMSE4);

tspan4 = linspace(t0_4, 36, 600);
f4 = @(t, V) a4 .* V .* (1 - V ./ b4);
[t_sol4, V_sol4] = ode45(f4, tspan4, V0_4);

figure();
plot(t_data4, V_data4, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold on;
plot(t_sol4, V_sol4, 'r-', 'LineWidth', 1.5);
hold off;

xlabel('Days');
ylabel('Tumor Volume (mm^3)');
title('Logistic BALB/Slc-nu/nu Mice');
xlim([5 37]);
ylim([0 2200]);
box on;

annotation('textbox', [0.3 0.01 0.4 0.08], 'String',sprintf('NMSE = %.7f\na = %.5f\nb = %.2f', NMSE4, a4, b4),'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);



%Function
function err = nmse_cost(params, t_data, V_data, t0, V0)
    a = params(1);
    b = params(2);
    if a <= 0 || b <= 0
        err = 1e6;
        return;
    end
    f = @(t, V) a .* V .* (1 - V ./ b);
    sol = ode45(f, [t0, t_data(end)], V0);
    V_model = deval(sol, t_data);
    err = sum((V_data - V_model).^2) / sum(V_data.^2);
end
