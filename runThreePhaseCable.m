clear; close all;
%% Parameters
addpath("Functions\");
addpath("Geometries\");
options.dmin = 0.3e-3;
options.gamma = 2;
options.N = 200;
options.linear = 0;

%% Three phase cable
p.x0 = 0;
p.y0 = 0;
p.r_w = 2e-3;                   % Wire radius
p.d_ins = 0.2e-3;               % Insulation thickness
p.d_w = 3e-3;                   % Distance between center of wire and cable center
p.L = 20e-3;                    % Box size

%% Parameters
p.theta = pi;
p.V0 = 100;                       % V
p.Rph = 2.5;                      % Ohm

p.Vu = p.V0*cos(p.theta);
p.Vv = p.V0*cos(p.theta-2*pi/3);
p.Vw = p.V0*cos(p.theta+2*pi/3);

p.I0 = p.V0/p.Rph;
p.Iu = p.I0*cos(p.theta);
p.Iv = p.I0*cos(p.theta-2*pi/3);
p.Iw = p.I0*cos(p.theta+2*pi/3);

Geom = returnThreePhaseCable(options,p);
[Geom] = linkNghbrBndrs(Geom);
Geom = returnNodes(Geom);

for nn = 1:length(Geom)
    [Geom(nn).G,Geom(nn).H] = calculateMatrixGHMultiBoundary(Geom(nn));
    [Geom(nn).G,Geom(nn).H,Geom(nn).bndr_val,Geom(nn).bndr_idx] = returnMultiDomainBndrVal(Geom(nn),Geom(nn).G,Geom(nn).H);
end
[A,B] = constructMultiDomain(Geom);
[Geom] = returnSolutionsMultiDomain(A,B,Geom);
%% Solve boundary values
figure
hold on
for kk = 1:length(Geom)
    for nn = 1:length(Geom(kk).boundary)
        I.idx0 = Geom(kk).boundary(nn).bndr_idx == 0;
        I.idx1 = Geom(kk).boundary(nn).bndr_idx == 1;
        I.idx2 = Geom(kk).boundary(nn).bndr_idx == 2;
        I.idx3 = Geom(kk).boundary(nn).bndr_idx == 3;

        plot(Geom(kk).boundary(nn).Xm(I.idx0),Geom(kk).boundary(nn).Ym(I.idx0),".k","MarkerSize",8)
        plot(Geom(kk).boundary(nn).Xm(I.idx1),Geom(kk).boundary(nn).Ym(I.idx1),".r","MarkerSize",8)
        plot(Geom(kk).boundary(nn).Xm(I.idx2),Geom(kk).boundary(nn).Ym(I.idx2),".g","MarkerSize",8)
        plot(Geom(kk).boundary(nn).Xm(I.idx3),Geom(kk).boundary(nn).Ym(I.idx3),".m","MarkerSize",8)
        addArrowLine(Geom(kk).boundary(nn).xgeom,Geom(kk).boundary(nn).ygeom)    
    end
end
grid on 
xlabel("x-axis [a.u.]","Interpreter","latex","FontSize",12)
ylabel("y-axis [a.u.]","Interpreter","latex","FontSize",12)

%% Solve domain values
Plt.x = p.x0 + 0.001e-3: 0.05e-3: p.x0 + p.L - 0.001e-3;
Plt.y = p.y0 + 0.001e-3: 0.05e-3: p.y0 + p.L - 0.001e-3;
[Plt.X, Plt.Y] = meshgrid(Plt.x, Plt.y); 

Plt.x01 = p.x0 + p.L/2; 
Plt.y01 = p.y0 + p.L/2 + p.d_w;
Plt.Idx1 = (Plt.X-Plt.x01).^2 + (Plt.Y-Plt.y01).^2 > p.r_w^2 & (Plt.X-Plt.x01).^2 + (Plt.Y-Plt.y01).^2 < (p.r_w + p.d_ins)^2;
Plt.Idx11 = (Plt.X-Plt.x01).^2 + (Plt.Y-Plt.y01).^2 <= p.r_w^2;

Plt.x02 = p.x0 + p.L/2 + cos(-pi/6)*p.d_w; 
Plt.y02 = p.y0 + p.L/2 + sin(-pi/6)*p.d_w;
Plt.Idx2 = (Plt.X-Plt.x02).^2 + (Plt.Y-Plt.y02).^2 > p.r_w^2 & (Plt.X-Plt.x02).^2 + (Plt.Y-Plt.y02).^2 < (p.r_w + p.d_ins)^2;
Plt.Idx22 = (Plt.X-Plt.x02).^2 + (Plt.Y-Plt.y02).^2 <= p.r_w^2;

Plt.x03 = p.x0 + p.L/2 + cos(-5*pi/6)*p.d_w; 
Plt.y03 = p.y0 + p.L/2 + sin(-5*pi/6)*p.d_w;
Plt.Idx3 = (Plt.X-Plt.x03).^2 + (Plt.Y-Plt.y03).^2 > p.r_w^2 & (Plt.X-Plt.x03).^2 + (Plt.Y-Plt.y03).^2 < (p.r_w + p.d_ins)^2;
Plt.Idx33 = (Plt.X-Plt.x03).^2 + (Plt.Y-Plt.y03).^2 <= p.r_w^2;

Plt.Idx4 = ~Plt.Idx1 & ~Plt.Idx11 & ~Plt.Idx2 & ~Plt.Idx22 & ~Plt.Idx3 & ~Plt.Idx33;

Uplot = zeros(size(Plt.X));


Uplot(Plt.Idx1) = calculateDomainValuesMulti(Plt.X(Plt.Idx1),Plt.Y(Plt.Idx1),Geom(1));
Uplot(Plt.Idx2) = calculateDomainValuesMulti(Plt.X(Plt.Idx2),Plt.Y(Plt.Idx2),Geom(2));
Uplot(Plt.Idx3) = calculateDomainValuesMulti(Plt.X(Plt.Idx3),Plt.Y(Plt.Idx3),Geom(3));
Uplot(Plt.Idx4) = calculateDomainValuesMulti(Plt.X(Plt.Idx4),Plt.Y(Plt.Idx4),Geom(4));

figure
hold on 

contourf(Plt.X,Plt.Y,Uplot,"FaceAlpha",0.5)
% quiver(Xplot(nn).val(1:dplot:end,1:dplot:end),Yplot(nn).val(1:dplot:end,1:dplot:end), ...
    % uplot_dx(nn).val(1:dplot:end,1:dplot:end),uplot_dy(nn).val(1:dplot:end,1:dplot:end),"LineWidth",1,"Color","k","AutoScale","on")
plot(Plt.X(Plt.Idx4),Plt.Y(Plt.Idx4),".k")

colormap("spring")
c = colorbar;
c.Label.String = 'Magnitude [a.u]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";

grid on
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("y-axis [m]","Interpreter","latex","FontSize",12)



% Check indices
% figure
% hold on 
% plot(Plt.X(Plt.Idx1),Plt.Y(Plt.Idx1),".r")
% plot(Plt.X(Plt.Idx2),Plt.Y(Plt.Idx2),".r")
% plot(Plt.X(Plt.Idx3),Plt.Y(Plt.Idx3),".r")
% plot(Plt.X(Plt.Idx4),Plt.Y(Plt.Idx4),".k")






