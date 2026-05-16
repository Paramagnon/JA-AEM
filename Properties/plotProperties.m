clear; close all;
%% Temp range
T = 4:0.1:300;                  % [K]
%% ---- Heat capacity
figure
loglog(T,Cp_alu(T),"LineWidth",2)
hold on 
loglog(T,Cp_copper(T),"LineWidth",2)
loglog(T,Cp_steel(T),"LineWidth",2)
legend("Aluminium","Copper","Stainless Steel","Location","best")
xlabel("Temperature [K]","Interpreter","latex","FontSize",12)
ylabel("Heat capacity [J/(kg.K)]","Interpreter","latex","FontSize",12)
grid on 

%% ---- Thermal conductivity
figure
loglog(T,kappa_alu(T,3003),"LineWidth",2)
hold on 
loglog(T,kappa_copper(T,50),"LineWidth",2)
loglog(T,kappa_steel(T),"LineWidth",2)
legend("Aluminium","Copper RRR 50","Stainless Steel","Location","best")
xlabel("Temperature [K]","Interpreter","latex","FontSize",12)
ylabel("Thermal conductivity [W/(m.K)]","Interpreter","latex","FontSize",12)
grid on 

%% ---- Electrical conductivity
figure
loglog(T,sigma_alu(T),"LineWidth",2)
hold on 
loglog(T,sigma_copper(T),"LineWidth",2)
loglog(T,sigma_steel(T),"LineWidth",2)
legend("Aluminium","Copper","Stainless Steel","Location","best")
xlabel("Temperature [K]","Interpreter","latex","FontSize",12)
ylabel("Electrical conductivity [S/m]","Interpreter","latex","FontSize",12)
grid on 

%% JcHTS
addpath("C:\Users\20182088\Documents\Graduation-Project\Matlab\AEM\Properties")
Tlow = 0:5:80; % K
Blow = 12:-0.5:0; % T
T = 0:0.1:80; % K
B = 12:-0.05:0; % T

JcHTS = zeros(length(T),length(B));
for tt = 1:length(T)
    for bb = 1:length(B)
        JcHTS(tt,bb) = returnJcHTS(T(tt),B(bb));
    end
end
[Tgrid,Bgrid] = meshgrid(B,T);

JcHTSlow = zeros(length(Tlow),length(Blow));
for tt = 1:length(Tlow)
    for bb = 1:length(Blow)
        JcHTSlow(tt,bb) = returnJcHTS(Tlow(tt),Blow(bb));
    end
end
[Tgridlow,Bgridlow] = meshgrid(Blow,Tlow);

addpath("C:\Users\20182088\Documents\Graduation-Project\Matlab\Functions")
figure(2)
surf(Tgrid,Bgrid,JcHTS*1e-9,"FaceColor","interp","EdgeColor","none")
hold on 
surf(Tgridlow,Bgridlow,JcHTSlow*1e-9,"FaceColor","interp","FaceColor","none","EdgeAlpha",0.5)
hold off
colorMapAltairNonLinear
%%
% colorMapThermal
ylabel("$T$ [K]","Interpreter","latex","FontSize",18)
xlabel("$B_\perp$ [T]","Interpreter","latex","FontSize",18)
zlabel("$J_{C}$ [kA/mm$^2$]","Interpreter","latex","FontSize",18)
ylim([T(1),T(end)])
xlim([B(end),B(1)])
zlim([0,returnJcHTS(T(1),0)*1e-9])
view(120,30)

%%
Ec = 1e-4;
n = 30;
Jc = 1.149e10;
J = linspace(0,1.4,100)*Jc;

EJ = @(j) sign(j).*((Ec*(abs(j)/Jc).^n)); 
Eplot = EJ(J);

figure(2)
plot(J/Jc,Eplot,"LineWidth",2)
grid on 
xlabel("Normalized current density J/Jc [-]","FontSize",12,"Interpreter","latex")
ylabel("Electric field [V/m]","FontSize",12,"Interpreter","latex")

%% Mover currents
tau_m = 0.03;
tau = 4*tau_m;
f = 20;
Nt = 300;
t = linspace(0,1/f,Nt);
dx = tau*sin(2*pi*f*t-pi/2)+tau;

figure(1)
% subplot(2,1,1)
yyaxis left
set(gca,'YColor','k');
plot(t*1000,dx,"LineWidth",1.5,"DisplayName","Position");
ylabel("Position [m]","FontSize",12,"Interpreter","latex")
yyaxis right
set(gca,'YColor','k');
plot(t*1000,gradient(dx,t),"LineWidth",1.5,"DisplayName","Velocity");
hold on
y = zeros(1,Nt);
y(1:Nt/2) = ones(1,Nt/2)*mean(gradient(dx(1:Nt/2),t(1:Nt/2)));
y(Nt/2+1:Nt) = ones(1,Nt/2)*mean(gradient(dx(Nt/2+1:Nt),t(Nt/2+1:Nt)));
% plot(t*1000,y,"LineWidth",1.5,"DisplayName",strcat("Average",string(newline),"velocity"))
plot(t*1000,y,"LineWidth",1.5,"DisplayName",strcat("Average"," velocity"))

hold off
grid on 
xlabel("Time [ms]","FontSize",12,"Interpreter","latex")
ylabel("Velocity [m/s]","FontSize",12,"Interpreter","latex")
set(gca,"Position",[0.1728 0.1680 0.6761 0.6023],"Units","normalized")
legend("NumColumns",3,"Units","normalized","Position",[0.08    0.8417    0.85    0.0793],"IconColumnWidth",20,"Interpreter","latex","FontSize",9)

Jref = 1e7;
Ja = Jref*cos(pi*dx/tau_m);
Jb = Jref*cos(pi*dx/tau_m-2*pi/3);
Jc = Jref*cos(pi*dx/tau_m+2*pi/3);

figure(2)
% subplot(2,1,2)
plot(t*1000,Ja,"LineWidth",1.5,"DisplayName","$J_a$")
hold on
plot(t*1000,Jb,"LineWidth",1.5,"DisplayName","$J_b$")
plot(t*1000,Jc,"LineWidth",1.5,"DisplayName","$J_c$")
hold off
grid on 
xlabel("Time [ms]","FontSize",12,"Interpreter","latex")
ylabel("Current density [A/m$^2$]","FontSize",12,"Interpreter","latex")
legend("Location","northoutside","NumColumns",3)
ax = gca();
ysecondarylabel(ax,'')
text(ax,'String','x10^{7}','Units','normalized','Position',[0 1.03],'HorizontalAlignment','right','VerticalAlignment','bottom')
ylim([-1.1e7 1.1e7])
set(gca,"Position",[0.1701 0.1680 0.8072 0.6023],"Units","normalized")
legend("Position",[0.2668    0.8417    0.6137    0.0793],"Units","normalized","Interpreter","latex","FontSize",9)


