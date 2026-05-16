clear; close all;
%% Parameters
addpath("Functions");
addpath("Functions\BEM\");
addpath("Functions\Mesh\");
addpath("Geometries\");
options.dmin = 0.5;
options.gamma = 1;
options.N = 2000;
options.linear = 1;
options.rmin = 1e-3;
options.dL = 100000;
%% Coldplate 
par.x0 = 0;
par.y0 = 0;
par.L = 200e-3;             % Length in [m]
par.H = 50e-3;              % Height in [m]
par.r = 5e-3;               % Radius in [m]

par.Npipe = 8;              % Number of pipes
par.T1 = 300;
par.T2 = 200;
%% Create Geom
[Geom] = returnColdPlate(options,par);
% [Geom] = returnTheTank(options,par);
figure
plotMesh(Geom)
domain = insertDomainPoints(Geom(1),options);

grid on 
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
axis square
%%
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
%% Test
tic 
x_unknown = A\B;
% [Geom] = linkNghbrBndrs(Geom);
% Geom = returnNodes(Geom);
% [Geom] = addDomainSource(Geom);
% [A,B] = constructMultiDomain(Geom);
% for nn = 1:length(Geom)
%     [Geom(nn).G,Geom(nn).H] = calculateMatrixGHMultiBoundary(Geom(nn));
%     [Geom(nn).G,Geom(nn).H,Geom(nn).bndr_val,Geom(nn).bndr_idx] = returnMultiDomainBndrVal(Geom(nn),Geom(nn).G,Geom(nn).H);
% end
% [Geom] = returnSolutionsMultiDomain(A,B,Geom);
toc
%% Domain solutions
x1 = linspace(par.x0+0.1e-3,par.L-0.1e-3,500); 
y1 = linspace(par.y0+0.1e-3,par.H-0.1e-3,300);
[Xplot,Yplot] = meshgrid(x1,y1);

par.dL = par.L/(par.Npipe+1);
idx = ones(size(Xplot));

for nn = 1:par.Npipe
    idx = idx & ((Xplot - par.x0 - nn*par.dL).^2 + (Yplot - par.y0 - par.H/2).^2 >= par.r^2);   
end

% figure
% plot(Xplot(idx),Yplot(idx),'.')

uplot = NaN(size(Xplot));

[uplot(idx)] = calculateDomainValuesMulti(Xplot(idx),Yplot(idx),Geom(1));

figure
hold on 
contourf(Xplot,Yplot,uplot,"FaceAlpha",0.7)

colormap("jet")
c = colorbar;
c.Label.String = 'Magnitude [K]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";
clim([0 par.T1+40])

grid on
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("y-axis [m]","Interpreter","latex","FontSize",12)