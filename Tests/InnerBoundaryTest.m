clear; close all;
%% Parameters
addpath("Functions\");
addpath("Functions\BEM\");
addpath("Functions\Mesh\");
addpath("Geometries\");
options.dmin = 0.3;
options.gamma = 1;
options.N = 20;
options.linear = 1;
%% Inner boundaries
close all
Geom = returnSquaresInSquares(options);

figure
hold on 
for kk = 1:length(Geom)
    for nn = 1:length(Geom(kk).boundary)
        plot(Geom(kk).boundary(nn).Xm,Geom(kk).boundary(nn).Ym,".-r","MarkerSize",8)
        addArrowLine(Geom(kk).boundary(nn).xgeom,Geom(kk).boundary(nn).ygeom)    
    end
end
xlabel("x-axis [a.u.]","Interpreter","latex","FontSize",12)
ylabel("y-axis [a.u.]","Interpreter","latex","FontSize",12)
grid on 
xlim([-0.5 5.5])
ylim([-0.5 5.5])
tic
%%
Geom = returnNodes(Geom);
for nn = 1:length(Geom)
    [Geom(nn).G,Geom(nn).H] = calculateMatrixGHMultiBoundary(Geom(nn));
    [Geom(nn).G,Geom(nn).H,Geom(nn).bndr_val,Geom(nn).bndr_idx] = returnMultiDomainBndrVal(Geom(nn),Geom(nn).G,Geom(nn).H);
end
[Geom(1).u_sol, Geom(1).un_sol] = returnSolutionsu(Geom(1).H,Geom(1).G,Geom(1).bndr_val,Geom(1).bndr_idx);
toc
%% Calculate domain values
% Square case
x = 0.001:0.03:4.999;
y = 0.001:0.03:4.999;
[Xplot,Yplot] = meshgrid(x,y);
[uplot] = calculateDomainValuesMulti(Xplot,Yplot,Geom);
[uplot_dx, uplot_dy] = calculateDomainValuesDerivativeMulti(Xplot,Yplot,Geom);

dplot = 10;
dumax = 200;
uplot_dx(abs(uplot_dx)>dumax) = NaN;
uplot_dy(abs(uplot_dy)>dumax) = NaN;

figure
contourf(Xplot,Yplot,uplot,"FaceAlpha",0.5)
colormap("summer")
c = colorbar;
c.Label.String = 'Magnitude [a.u]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";
hold on
quiver(Xplot(1:dplot:end,1:dplot:end),Yplot(1:dplot:end,1:dplot:end),uplot_dx(1:dplot:end,1:dplot:end),uplot_dy(1:dplot:end,1:dplot:end),1,"LineWidth",1,"Color","k")
grid on
for nn = 2:4
    P1 = polyshape([Geom(1).boundary(nn).xgeom(1:end-1)],[Geom(1).boundary(nn).ygeom(1:end-1)]);
    plot(P1,"FaceAlpha",1,"FaceColor","k")
end
circles(Geom(1).boundary(5).p.x0,Geom(1).boundary(5).p.y0,Geom(1).boundary(5).p.r,"facecolor","black")

xlabel("x-axis [a.u.]","Interpreter","latex","FontSize",12)
ylabel("y-axis [a.u.]","Interpreter","latex","FontSize",12)
xlim([-0.5 5.5])
ylim([-0.5 5.5])