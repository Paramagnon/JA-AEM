clear; close all;
%% Parameters
addpath("Functions");
addpath("Functions\BEM\");
addpath("Functions\Mesh\");
addpath("Geometries\");
options.dmin = 0.5;                 % Min distance from points
options.gamma = 1;                  % Set exponential mesh density
options.N = 2000;                   % Specifies a linear mesh density
options.linear = 1;                 % Toggle between mesh modes

%% Magnet test
Geom = returnMagnetTest(options);
figure
plotMesh(Geom)
grid on 
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("y-axis [m]","Interpreter","latex","FontSize",12)

[Geom] = linkNghbrBndrs(Geom);
Geom = returnNodes(Geom);
[Geom] = addDomainSource(Geom);

tic
for nn = 1:length(Geom)
    [Geom(nn).G,Geom(nn).H] = calculateMatrixGHMultiBoundary(Geom(nn));
    [Geom(nn).G,Geom(nn).H,Geom(nn).bndr_val,Geom(nn).bndr_idx] = returnMultiDomainBndrVal(Geom(nn),Geom(nn).G,Geom(nn).H);
end
[A,B] = constructMultiDomain(Geom);
[Geom] = returnSolutionsMultiDomain(A,B,Geom);
toc
%% Domain solutions
% Square 1
x1 = [0.001:0.5:119.999]*1e-3; 
y1 = [10.001:0.5:34.999]*1e-3;
[Xplot(1).val,Yplot(1).val] = meshgrid(x1,y1);

% x2 = 5.001:0.03:9.999;
% y2 = 0.001:0.03:4.999;
% [Xplot(2).val,Yplot(2).val] = meshgrid(x2,y2);
% 
% x3 = 0.001:0.03:4.999;
% y3 = 5.001:0.03:9.999;
% [Xplot(3).val,Yplot(3).val] = meshgrid(x3,y3);

dplot = 10;
dumaxx = 0.1;
dumaxy = 1;

for nn = 1
    [uplot(nn).val] = calculateDomainValuesMulti(Xplot(nn).val,Yplot(nn).val,Geom(5));
    [uplot_dx(nn).val, uplot_dy(nn).val] = calculateDomainValuesDerivativeMulti(Xplot(nn).val,Yplot(nn).val,Geom(5));
    uplot_dx(nn).val(abs(uplot_dx(nn).val)>dumaxx) = NaN;
    uplot_dy(nn).val(abs(uplot_dy(nn).val)>dumaxy) = NaN;
end

figure
hold on 
for nn = 1
    contourf(Xplot(nn).val,Yplot(nn).val,uplot(nn).val,"FaceAlpha",0.5)
    quiver(Xplot(nn).val(1:dplot:end,1:dplot:end),Yplot(nn).val(1:dplot:end,1:dplot:end), ...
        uplot_dy(nn).val(1:dplot:end,1:dplot:end),-uplot_dx(nn).val(1:dplot:end,1:dplot:end),"LineWidth",1,"Color","k","AutoScale","on")
end
colormap("summer")
c = colorbar;
c.Label.String = 'Magnitude [a.u]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";

grid on
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("y-axis [m]","Interpreter","latex","FontSize",12)

%% Plot line potential derivatives
xline = [0.1:0.5:119.9]*1e-3;
% xline = [0.001:0.01:47.999]*1e-3;
mu0 = 4*pi*1e-7';

yline = ones(size(xline))*25e-3;
[Hx, Hy] = calculateDomainValuesDerivativeMulti(xline,yline,Geom(5));
%--- For magnetic potential
% Bx = mu0*Hx;
% By = mu0*Hy;
%--- For magnetic vector potential
Bx = Hy;
By = -Hx;

figure
plot(xline,Bx,"LineWidth",2)
hold on 
plot(xline,By,"LineWidth",2)
grid on 
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("B field [T]","Interpreter","latex","FontSize",12)
legend("Bx","By","Location","best")

%% Plot line potential
xline = [0.001:0.5:119.999]*1e-3;
yline = [];
yline(1).val = ones(size(xline))*12e-3;
yline(2).val = ones(size(xline))*14e-3;
yline(3).val = ones(size(xline))*18e-3;

figure
hold on
for nn = 1:3
    [uline(nn).val] = calculateDomainValuesMulti(xline,yline(nn).val,Geom(5));
    plot(xline,uline(nn).val,"LineWidth",2)   
end
legend("h = 8e-3","h = 10e-3","h = 14e-3","Location","best")

grid on 
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("Potential u [Wb/m]","Interpreter","latex","FontSize",12)



