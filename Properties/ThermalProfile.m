% --- Parameters ---
L1 = 200e-3;        % Length of segment 1
L2 = 200e-3;        % Length of segment 2
dwall = 2e-3;

k1 = kappa_steel(290,1);       % Conductivity segment 1
k2 = kappa_steel(290,1);        % Conductivity segment 2

q0 = 0;      % Heat flux at x = 0
Twall = 293;
qrad = 0.2*5.67e-8*Twall^4;
qeddy = 0.1/(L1*dwall); 
qgen = qeddy + qrad;      % Heat generation ONLY in segment 1
qgen = 3e6;
TR = 77;      % Right boundary temperature

% --- Discretization ---
N = 200;
x = linspace(0, L1 + L2, N);

% --- Constants ---
C1 = -q0 / k1;

C3 = (-qgen * L1 - q0) / k2;

C4 = TR - C3 * (L1 + L2);

C2 = C3 * L1 + C4 ...
     + (qgen/(2*k1)) * L1^2 ...
     + (q0/k1) * L1;

% --- Temperature profile ---
T = zeros(size(x));

for i = 1:length(x)
    if x(i) <= L1
        % Segment 1 (with generation)
        T(i) = -(qgen/(2*k1))*x(i)^2 ...
               - (q0/k1)*x(i) ...
               + C2;
    else
        % Segment 2 (no generation)
        T(i) = C3 * x(i) + C4;
    end
end

% --- Plot ---
figure(1);
subplot(2,1,1)
plot(x, T, 'LineWidth', 2);
hold on
xline(L1, '--r', 'Interface');
hold off
grid on
xlabel('x',"FontSize",12,"Interpreter","latex");
ylabel('Temperature',"FontSize",12,"Interpreter","latex");
legend('Temperature', 'Interface');

% Heat flux at right boundary
q_out = qgen * L1 + q0;

% If cross-sectional area is known
A = 1.0; % example
Q_out = A * q_out;

fprintf('Heat flux at right boundary: %.2f W/m^2\n', q_out);
fprintf('Total heat removal required: %.2f W\n', Q_out);

emiss = 0.2;
T1 = T(x<L1);
qradHe = 5.67e-8*emiss*(T1.^4-4^4);
qradH2 = 5.67e-8*emiss*(T1.^4-20^4);
qradN2 = 5.67e-8*emiss*(T1.^4-77^4);

qshld = 5.67e-8*emiss*(293^4-T1.^4);

subplot(2,1,2)
plot(x(x<L1),qradHe,"LineWidth",2)
hold on 
plot(x(x<L1),qradH2,"LineWidth",2)
plot(x(x<L1),qradN2,"LineWidth",2)
% plot(x(x<L1),qshld,"LineWidth",2)
hold off
grid on 
xlabel('x',"FontSize",12,"Interpreter","latex");
ylabel('Radiation heat flux [W/m$^2$]',"FontSize",12,"Interpreter","latex");
legend("LHe","LH2","LN2")
xlim([x(1), x(end)])
%% Skin depth
T0 = 100;
f = 100;
dlta = 1000*sqrt(1/(sigma_copper(T0)*pi*f*4*pi*1e-7))

% dlam = 1e-3;
% Ba = 0.4;
% Pe = pi^2*dlam^2*f^2*Ba^2*sigma_copper(T0)/6
