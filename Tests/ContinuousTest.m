clear; close all;
%% Parameters
addpath("Functions");
addpath("Functions\BEM\");
addpath("Functions\Mesh\");
addpath("AEM\TimeSolver\");
addpath("AEM\Matrices\");
addpath("AEM\Mesh\");
addpath("Geometries\");
options.dmin = 1e-3;
options.gamma = 1;
options.N = 2000;
options.linear = 1;
options.dL = 400000;
options.rmin = 0.001;

%% Continuous test
[Geom] = returnContinuousTest(options);
[Geom] = linkNghbrBndrs(Geom);
Geom = returnNodes(Geom);

SimStep = runTime(Geom);
alpha = SimStep.K\SimStep.p;


%% Figures
src = @(t) -2e6*sin(2*pi*10*t);

L1 = Geom(1).Ldmn;
SimStep.f = zeros(length(SimStep.S00),1);
SimStep.f(L1+1:end) = 1;

Cidx = zeros(size(SimStep.Uxx));
Cidx(1:L1,1:L1) = eye(L1);
C = Cidx*SimStep.S00;
K = (SimStep.Uxx + SimStep.Uyy);

tic 
dt = 0.001;
t = 0:dt:0.1;


alpha_mat = zeros(length(SimStep.f),length(t));
dadt = zeros(length(SimStep.f),length(t));

for tt = 1:length(t)-1
    f =  SimStep.f*src(t(tt+1));

    alpha_mat(:,tt+1) = (C+dt*K)\(dt*f + C*alpha_mat(:,tt));
    % alpha_mat(:,tt) = (SimStep.Uxx+SimStep.Uyy)\f;

    dadt(:,tt+1) = (SimStep.S00*alpha_mat(:,tt+1)-SimStep.S00*alpha_mat(:,tt))/dt;
end

toc 



%% ---------- Magnitude
figure
colormap("jet")
for n = 1:length(t)       
    plotFrameA(alpha_mat(:,n), SimStep, Geom)
    c = colorbar;
    c.Label.String = 'Magnitude [Wb/m]';
    c.Label.FontSize = 12;
    c.Label.Interpreter = "latex";
    
    clim([-25 25])

    grid on
    xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
    ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
    axis square

    exportgraphics(gcf,'testTime.gif','Append',true);
end

%% ---------- Power loss
figure
colormap("jet")
for n = 1:length(t)       
    plotFramedAdt(Cidx*dadt(:,n).^2, Geom)
    c = colorbar;
    c.Label.String = 'Power [W]';
    c.Label.FontSize = 12;
    c.Label.Interpreter = "latex";
    
    clim([0 8e5])

    grid on
    xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
    ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
    axis square

    exportgraphics(gcf,'testTimePower.gif','Append',true);
end

%%

% --------------- Plot mesh
% figure
% plotMesh(Geom)
% axis square
% grid on 
% xlabel("x-axis [a.u.]","Interpreter","latex","FontSize",12)
% ylabel("y-axis [a.u.]","Interpreter","latex","FontSize",12)

% --------------- plot contour
u = SimStep.S00*alpha+SimStep.z00;
ux = SimStep.Sx0*alpha+SimStep.zx0;
uy = SimStep.Sy0*alpha+SimStep.zy0;

L = 0;
for kk = 1:length(Geom)    
    L = L + Geom(kk).Ldmn;    
end
Npart = length(Geom(1).boundary(1).bndr_idx);
x = zeros(L,1);
y = zeros(L,1);
L = L + Npart;
    
xp = zeros(L,1);
yp = zeros(L,1);

xp(1:Npart) = Geom(1).boundary(1).Xn;
yp(1:Npart) = Geom(1).boundary(1).Yn;

up = zeros(L,1);
uq = SimStep.Z*alpha  + SimStep.e;

up(1:Npart) = uq(1:Npart);
up(Npart+1:L) = u;

Loffset = 0;
for kk = 1:length(Geom)        
    Lpart = Geom(kk).Ldmn;
    x(Loffset+1:Loffset+Lpart) = Geom(kk).domain.Xn;
    y(Loffset+1:Loffset+Lpart) = Geom(kk).domain.Yn;
    Loffset = Loffset + Lpart;
end
xp(Npart+1:L) = x;
yp(Npart+1:L) = y;

xv = linspace(min(xp), max(xp), 200);
yv = linspace(min(yp), max(yp), 200);

[xq, yq] = ndgrid(xv, yv);
maggrid = griddata(xp(:), yp(:), up(:), xq, yq);

xv2 = linspace(min(x), max(x), 10);
yv2 = linspace(min(y), max(y), 10);

[xq2, yq2] = ndgrid(xv2, yv2);


maggridx = griddata(x(:), y(:), ux(:), xq2, yq2);
maggridy = griddata(x(:), y(:), uy(:), xq2, yq2);

figure
contourf(xq, yq, maggrid, 20);
hold on
quiver(xq2,yq2,-maggridy,maggridx,"LineWidth",1,"Color","k","AutoScale","on")
colormap("summer")
c = colorbar;
c.Label.String = 'Magnitude [Wb/m]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";

grid on
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
axis square

% --------------- plot u
xline = 0.015:0.00005:0.03;
yline = 0.015*ones(size(xline));

maggridline = griddata(xp(:), yp(:), up(:), xline, yline);

figure
plot(xline,maggridline,"LineWidth",2)
grid on 
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("Magnitude [Wb/m]","Interpreter","latex","FontSize",12)

%% Functions

function plotFrameA(alpha, SimStep, Geom)
    u = SimStep.S00*alpha+SimStep.z00;
    L = 0;
    for kk = 1:length(Geom)    
        L = L + Geom(kk).Ldmn;    
    end
    Npart = length(Geom(1).boundary(1).bndr_idx);
    x = zeros(L,1);
    y = zeros(L,1);
    L = L + Npart;
        
    xp = zeros(L,1);
    yp = zeros(L,1);
    
    xp(1:Npart) = Geom(1).boundary(1).Xn;
    yp(1:Npart) = Geom(1).boundary(1).Yn;
    
    up = zeros(L,1);
    uq = SimStep.Z*alpha  + SimStep.e;
    
    up(1:Npart) = uq(1:Npart);
    up(Npart+1:L) = u;
    
    Loffset = 0;
    for kk = 1:length(Geom)        
        Lpart = Geom(kk).Ldmn;
        x(Loffset+1:Loffset+Lpart) = Geom(kk).domain.Xn;
        y(Loffset+1:Loffset+Lpart) = Geom(kk).domain.Yn;
        Loffset = Loffset + Lpart;
    end
    xp(Npart+1:L) = x;
    yp(Npart+1:L) = y;
    
    xv = linspace(min(xp), max(xp), 200);
    yv = linspace(min(yp), max(yp), 200);
    
    [xq, yq] = ndgrid(xv, yv);
    maggrid = griddata(xp(:), yp(:), up(:), xq, yq);
    
    contourf(xq, yq, maggrid, 20);
    
end

function plotFramedAdt(uy, Geom)
    L = 0;
    for kk = 1:length(Geom)    
        L = L + Geom(kk).Ldmn;    
    end    
    x = zeros(L,1);
    y = zeros(L,1);
   
    Loffset = 0;
    for kk = 1:length(Geom)        
        Lpart = Geom(kk).Ldmn;
        x(Loffset+1:Loffset+Lpart) = Geom(kk).domain.Xn;
        y(Loffset+1:Loffset+Lpart) = Geom(kk).domain.Yn;
        Loffset = Loffset + Lpart;
    end 

    xv2 = linspace(min(x), max(x), 200);
    yv2 = linspace(min(y), max(y), 200);

    [xq2, yq2] = ndgrid(xv2, yv2);

    maggridy = griddata(x(:), y(:), uy(:), xq2, yq2);
    
    contourf(xq2, yq2, maggridy, 20);

end

