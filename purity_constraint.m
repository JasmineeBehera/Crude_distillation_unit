function [c, ceq] = purity_constraint(R)
    Nstage = 8;
    feedTray = 4;
    alpha = 2.5;
    F = 100;
    zF = 0.5;
    D = 50;
    x0 = linspace(0.9,0.1,Nstage)';

    options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',5000,'MaxIterations',2000);
    x = fsolve(@(x) tray_mass_energy(x,Nstage,alpha,F,zF,R,D,feedTray), x0, options);
    x = max(min(x,0.999),0.001);

    x_top = x(1);

    c = 0.95 - x_top; % top purity ≥ 0.95
    ceq = [];
end
