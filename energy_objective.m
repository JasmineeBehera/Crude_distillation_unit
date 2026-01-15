function E_total = energy_objective(R)
    % Column parameters
    Nstage = 8;
    feedTray = 4;
    alpha = 2.5;
    F = 100;
    zF = 0.5;
    D = 50;

    % Initial guess
    x0 = linspace(0.9,0.1,Nstage)';

    % Solve tray balances
    options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',5000,'MaxIterations',2000);
    x = fsolve(@(x) tray_mass_energy(x,Nstage,alpha,F,zF,R,D,feedTray), x0, options);

    % Physical bounds
    x = max(min(x,0.999),0.001);

    % Energy calculation
    Hvap = 35000;  
    V = R + 1;     

    x_top = x(1);
    x_bottom = x(end);

    Q_reboiler = V * Hvap * x_bottom;
    Q_condenser = V * Hvap * x_top;

    E_total = Q_reboiler + Q_condenser;
end
