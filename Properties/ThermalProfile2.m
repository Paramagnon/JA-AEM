% ============================
% PARAMETERS
% ============================

L1 = 400e-3;              % Segment 1 length [m]
L2 = 300e-3;              % Segment 2 length [m]

k1 = 10.0;             % Conductivity segment 1 [W/mK]
k2 = 5.0;              % Conductivity segment 2 [W/mK]

q0 = 0;            % Heat flux at x=0 [W/m^2]
q_eddy = 500;          % Eddy volumetric heating [W/m^3]

eps = 0.2;             % Emissivity
sigma = 5.67e-8;       % Stefan-Boltzmann [W/m^2K^4]
t = 0.01;              % Thickness [m]
T_inf = 293;           % محیط (shield temp) [K]

TR = 80;               % Cold end wall temperature guess [K]

A_cross = 1e-4;        % Cross-sectional area [m^2]

% Cryogen properties (example: liquid nitrogen)
cp = 2000;             % J/kg-K
T_fluid_in = 77;       % K
dT_fluid_max = 5;      % Allowed ताप rise [K]

% ============================
% DISCRETIZATION
% ============================

N = 200;
x = linspace(0, L1+L2, N);
dx = x(2)-x(1);

T = ones(N,1)*T_inf;   % initial guess

% ============================
% NONLINEAR SOLVER (FDM)
% ============================

max_iter = 5000;
tol = 1e-6;

for iter = 1:max_iter
    T_old = T;
    
    for i = 2:N-1
        
        if x(i) <= L1
            % Segment 1: conduction + eddy + radiation
            
            qrad = eps * sigma * (T_inf^4 - T(i)^4);
            
            source = q_eddy + qrad;
            
            T(i) = 0.5 * (T(i+1) + T(i-1) + dx^2 * source / k1);
            
        else
            % Segment 2: pure conduction
            
            T(i) = 0.5 * (T(i+1) + T(i-1));
        end
    end
    
    % Left boundary (flux)
    T(1) = T(2) + q0 * dx / k1;
    
    % Right boundary (fixed temperature)
    T(end) = TR;
    
    % Convergence check
    if max(abs(T - T_old)) < tol
        fprintf('Converged in %d iterations\n', iter);
        break;
    end
end

% ============================
% HEAT FLUX AT COLD END
% ============================

q_out = -k2 * (T(end) - T(end-1)) / dx;   % W/m^2
Q_out = q_out * A_cross;                 % Total heat [W]

fprintf('Heat flux out: %.2f W/m^2\n', q_out);
fprintf('Total heat load: %.4f W\n', Q_out);

% ============================
% REQUIRED CONVECTION COEFFICIENT
% ============================

h_required = q_out / (T(end) - T_fluid_in);

fprintf('Required h: %.2f W/m^2K\n', h_required);

% ============================
% MASS FLOW RATE CALCULATION
% ============================

% Energy balance: Q = m_dot * cp * dT
m_dot = Q_out / (cp * dT_fluid_max);

fprintf('Required mass flow rate: %.6f kg/s\n', m_dot);

% ============================
% PLOT
% ============================

figure(2);
plot(x, T, 'LineWidth', 2);
xlabel('x [m]');
ylabel('Temperature [K]');
title('Temperature Profile with Radiation + Cooling');
grid on;
xline(L1, '--r', 'Interface');