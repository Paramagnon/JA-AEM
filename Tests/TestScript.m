clear; close all;
%% Parameters
Jc0 = 2.85e10;          % [A/m^2]
Ec = 1e-4;
mu0 = 4*pi*1e-7;
Aw = 2e-6*4e-3;         % 2 um * 4 mm
%--- Material properties
alpha = 1.5; 
beta = 1;
gamma = 5.02;
B0 = 240e-3;            % [T]
kappa = 4;
n = 30;
%% Test
J = logspace(1,5,500);

nablaAl = -(mu0*1-kappa*Aw)*J;
nablaAr = Ec*sign(J).*abs((J/Jc0).^n);
figure
semilogx(J,-nablaAl)
hold on 
semilogx(J,nablaAr)
grid on 
xlabel("Current density [A/m^2]")
ylabel("\nabla^2 A")


