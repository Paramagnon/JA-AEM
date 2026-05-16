%% Heat transfer coefficient
par.r_pipes = 1e-3;
options.cryoflow = 1;
%----- Heat transfer for different cryogens
gm.A = pi*par.r_pipes^2;                   % [m^2]
gm.u = 0.001*options.cryoflow/gm.A;           % [m/s]
gm.L = 2*par.r_pipes;                      % [m] Diameter of the pipes
% Tcryo = par.Tcryo;

options.cryotype = 1;
T1 = linspace(63.2,77.3,100);
h1 = zeros(size(T1));
for tt = 1:length(T1)
    h1(tt) = returnh(T1(tt),gm,options);
end

options.cryotype = 2;
T2 = linspace(14,20.3,100);
h2 = zeros(size(T2));
for tt = 1:length(T2)
    h2(tt) = returnh(T2(tt),gm,options);
end

options.cryotype = 3;
T3 = linspace(2.2,4.2,100);
h3 = zeros(size(T3));
for tt = 1:length(T3)
    h3(tt) = returnh(T3(tt),gm,options);
end

figure
hold on 
plot(T1,h1,"LineWidth",2)
plot(T2,h2,"LineWidth",2)
plot(T3,h3,"LineWidth",2)

grid on 
xlabel("Temperature [K]","FontSize",12,"Interpreter","latex")
ylabel("Heat transfer [W/(m$^2$.k)]","FontSize",12,"Interpreter","latex")
%% Properties
figure(2)
subplot(2,2,1)
hold on
plot(T1,kappa_lnitrogen(T1),"LineWidth",2)
plot(T2,kappa_lhydrogen(T2),"LineWidth",2)
plot(T3,kappa_lhelium(T3),"LineWidth",2)
grid on 
legend("Nitrogen","Hydrogen","Helium","Location","northwest")
xlabel("Temperature [K]","FontSize",12,"Interpreter","latex")
ylabel("Thermal conductivity [W/(m.k)]","FontSize",12,"Interpreter","latex")

subplot(2,2,2)
hold on
plot(T1,rho_lnitrogen(T1),"LineWidth",2)
plot(T2,rho_lhydrogen(T2),"LineWidth",2)
plot(T3,rho_lhelium(T3),"LineWidth",2)
grid on 
legend("Nitrogen","Hydrogen","Helium","Location","northwest")
xlabel("Temperature [K]","FontSize",12,"Interpreter","latex")
ylabel("Density [kg/m$^3$)]","FontSize",12,"Interpreter","latex")

subplot(2,2,3)
hold on
plot(T1,Cp_lnitrogen(T1),"LineWidth",2)
plot(T2,Cp_lhydrogen(T2),"LineWidth",2)
plot(T3,Cp_lhelium(T3),"LineWidth",2)
grid on 
legend("Nitrogen","Hydrogen","Helium","Location","northeast")
xlabel("Temperature [K]","FontSize",12,"Interpreter","latex")
ylabel("Specific heat [J/(kg.k)]","FontSize",12,"Interpreter","latex")

subplot(2,2,4)
hold on
plot(T1,muv_lnitrogen(T1),"LineWidth",2)
plot(T2,muv_lhydrogen(T2),"LineWidth",2)
plot(T3,muv_lhelium(T3),"LineWidth",2)
grid on 
legend("Nitrogen","Hydrogen","Helium","Location","northwest")
xlabel("Temperature [K]","FontSize",12,"Interpreter","latex")
ylabel("Dynamic viscosity [Pa.s]","FontSize",12,"Interpreter","latex")

%% Thermal properties
T = 4:0.1:80;

figure
% ---- Heat capacity
subplot(2,2,1)
semilogy(T,Cp_alu(T),"LineWidth",2)
hold on 
semilogy(T,Cp_copper(T),"LineWidth",2)
semilogy(T,Cp_steel(T),"LineWidth",2)
legend("Aluminium","Copper","Stainless Steel","Location","best")
xlabel("Temperature [K]","Interpreter","latex","FontSize",12)
ylabel("Heat capacity [J/(kg.K)]","Interpreter","latex","FontSize",12)
grid on 

% ---- Thermal conductivity
subplot(2,2,2)
semilogy(T,kappa_alu(T,3003),"LineWidth",2)
hold on 
semilogy(T,kappa_copper(T,50),"LineWidth",2)
semilogy(T,kappa_steel(T),"LineWidth",2)
legend("Aluminium","Copper RRR 50","Stainless Steel","Location","best")
xlabel("Temperature [K]","Interpreter","latex","FontSize",12)
ylabel("Thermal conductivity [W/(m.K)]","Interpreter","latex","FontSize",12)
grid on 

% ---- Electrical conductivity
subplot(2,2,3)
semilogy(T,sigma_alu(T),"LineWidth",2)
hold on 
semilogy(T,sigma_copper(T),"LineWidth",2)
semilogy(T,sigma_steel(T),"LineWidth",2)
legend("Aluminium","Copper","Stainless Steel","Location","best")
xlabel("Temperature [K]","Interpreter","latex","FontSize",12)
ylabel("Electrical conductivity [S/m]","Interpreter","latex","FontSize",12)
grid on 

% ---- Density
a_alu = kappa_alu(T,3003)./(rho_alu(T).*Cp_alu(T));
a_copper = kappa_copper(T,50)./(rho_copper(T).*Cp_copper(T));
a_steel = kappa_steel(T)./(rho_steel(T).*Cp_steel(T));

subplot(2,2,4)
semilogy(T,a_alu,"LineWidth",2)
hold on 
semilogy(T,a_copper,"LineWidth",2)
semilogy(T,a_steel,"LineWidth",2)
legend("Aluminium","Copper","Stainless Steel","Location","best")
xlabel("Temperature [K]","Interpreter","latex","FontSize",12)
ylabel("Thermal diffusivity [m$^2$/s]","Interpreter","latex","FontSize",12)
grid on 

%% Heat transfer multi plate
SB = 5.67e-8; % Boltzmann constant0
T1 = 293; % [k]
T2 = 80; % [k]
T3 = 4; % [K]

em1 = 0.2;
em2 = 0.2;
em3 = 0.3;

q1 = SB*(T1^4-T2^4)/(1/em1 + 1/em2 - 1)
q2 = SB*(T2^4-T3^4)/(1/em2 + 1/em3 - 1)

%% Functions
function h = returnh(Tcryo,gm,options)   
    switch options.cryotype
        case 1 % Liquid Nitrogen
            if Tcryo > 77.355 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] above boiling point'])
            elseif Tcryo < 63.151 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] below data range'])
            end
            Pr = Cp_lnitrogen(Tcryo)*muv_lnitrogen(Tcryo)/kappa_lnitrogen(Tcryo);
            Re = rho_lnitrogen(Tcryo)*gm.u*gm.L/muv_lnitrogen(Tcryo);
            Nu = 0.027*Re^0.8*Pr^0.4;
            h = Nu*kappa_lnitrogen(Tcryo)/gm.L;
  
        case 2 % Liquid Hydrogen           
            if Tcryo > 20.369 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] above boiling point'])
            elseif Tcryo < 13.957 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] below data range'])
            end
            Pr = Cp_lhydrogen(Tcryo)*muv_lhydrogen(Tcryo)/kappa_lhydrogen(Tcryo);
            Re = rho_lhydrogen(Tcryo)*gm.u*gm.L/muv_lhydrogen(Tcryo);
            Nu = 0.023*Re^0.8*Pr^0.4;
            h = Nu*kappa_lhydrogen(Tcryo)/gm.L;

        case 3 % Liquid Helium            
            if Tcryo > 4.2238 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] above boiling point'])
            elseif Tcryo < 2.177 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] below data range'])
            end
            Pr = Cp_lhelium(Tcryo)*muv_lhelium(Tcryo)/kappa_lhelium(Tcryo);
            Re = rho_lhelium(Tcryo)*gm.u*gm.L/muv_lhelium(Tcryo);
            Nu = 0.022*Re^0.8*Pr^0.4;
            h = Nu*kappa_lhelium(Tcryo)/gm.L;

        otherwise
            disp("Not a valid cryogen type");
    end
end