function Fbal = tray_mass_energy(x,N,alpha,F,zF,R,D,feedTray)

% Physical bounds (ESSENTIAL for numerical stability)
x = max(min(x,0.999),0.001);

Fbal = zeros(N,1);

% VLE
y = (alpha.*x)./(1 + (alpha-1).*x);

% Normalized internal flows
L = R;
V = R + 1;

for i = 1:N

    if i == 1
        % Total condenser
        Fbal(i) = V*y(i+1) - (L+1)*x(i);

    elseif i == feedTray
        % Feed tray
        Fbal(i) = ...
            L*x(i-1) + V*y(i+1) + F*zF ...
            - (L+F)*x(i) - V*y(i);

    elseif i == N
        % Reboiler
        Fbal(i) = (L+F)*x(i-1) - V*y(i);

    else
        % Internal trays
        Fbal(i) = ...
            L*x(i-1) + V*y(i+1) ...
            - L*x(i) - V*y(i);
    end
end
end
