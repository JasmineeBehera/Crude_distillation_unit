clc; clear;

%% ---- Column parameters ----
Nstage = 16;      % increased trays for feasible purity
feedTray = 7;

alpha = 2.5;      % relative volatility
F = 100;           % feed flow (mol/s)
zF = 0.5;          % feed composition (light key)
D = 50;            % distillate flow
R0 = 2;            % initial guess for reflux

%% ---- Initial guess for tray compositions ----
x0 = linspace(0.9,0.1,Nstage)';

options_fsolve = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',5000,'MaxIterations',2000);

%% ---- Solve tray balances for initial reflux ----
x = fsolve(@(x) tray_mass_energy(x,Nstage,alpha,F,zF,R0,D,feedTray), x0, options_fsolve);
x = max(min(x,0.999),0.001);

disp('Tray compositions (light key):')
disp(x)

%% ---- Initial Energy Calculation ----
Hvap = 35000;  
V = R0 + 1;  

x_top = x(1);
x_bottom = x(end);

Q_reboiler = V * Hvap * x_bottom;
Q_condenser = V * Hvap * x_top;
Total_energy = Q_reboiler + Q_condenser;

disp(['Total energy (J/mol) = ', num2str(Total_energy)])

%% ---- Optimization using fmincon ----
lb = 0.5; ub = 5;  % reflux bounds

options_opt = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
    'StepTolerance',1e-8, ...
    'ConstraintTolerance',1e-4, ...   % <-- slightly relaxed
    'OptimalityTolerance',1e-6, ...
    'MaxFunctionEvaluations',5000);


[R_opt, E_min] = fmincon(@(R) energy_objective(R), R0, [], [], [], [], lb, ub, @(R) purity_constraint(R), options_opt);

disp(['Optimal Reflux Ratio = ', num2str(R_opt)])
disp(['Minimum Total Energy = ', num2str(E_min)])

%% ---- Energy vs Reflux Ratio Plot ----
R_range = 0.5:0.2:5;
Energy = zeros(length(R_range),1);
Top_purity = zeros(length(R_range),1);

for k = 1:length(R_range)
    R = R_range(k);
    x = fsolve(@(x) tray_mass_energy(x,Nstage,alpha,F,zF,R,D,feedTray),x0,options_fsolve);
    x = max(min(x,0.999),0.001);

    V = R + 1;
    x_top = x(1);
    x_bottom = x(end);

    Energy(k) = V*Hvap*x_top + V*Hvap*x_bottom;
    Top_purity(k) = x_top;
end

figure
yyaxis left
plot(R_range,Energy,'-o','LineWidth',2)
ylabel('Total Energy (J/mol)')
yyaxis right
plot(R_range,Top_purity,'-s','LineWidth',2)
xlabel('Reflux Ratio')
ylabel('Top Product Purity')
grid on
legend('Total Energy','Top Purity')
title('Energy and Purity vs Reflux Ratio')

%% --- Local Functions ---
function E_total = energy_objective(R)
    Nstage = 12; feedTray = 5; alpha = 2.5; F = 100; zF = 0.5; D = 50;
    x0 = linspace(0.9,0.1,Nstage)';
    options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',5000,'MaxIterations',2000);
    x = fsolve(@(x) tray_mass_energy(x,Nstage,alpha,F,zF,R,D,feedTray), x0, options);
    x = max(min(x,0.999),0.001);
    Hvap = 35000; V = R+1;
    x_top = x(1); x_bottom = x(end);
    Q_reboiler = V*Hvap*x_bottom; Q_condenser = V*Hvap*x_top;
    E_total = Q_reboiler + Q_condenser;
end

function [c, ceq] = purity_constraint(R)
    Nstage = 12; feedTray = 5; alpha = 2.5; F = 100; zF = 0.5; D = 50;
    x0 = linspace(0.9,0.1,Nstage)';
    options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',5000,'MaxIterations',2000);
    x = fsolve(@(x) tray_mass_energy(x,Nstage,alpha,F,zF,R,D,feedTray), x0, options);
    x = max(min(x,0.999),0.001);
    x_top = x(1);
    c = 0.949 - x_top; % top purity ≥ 0.95
    ceq = [];
end