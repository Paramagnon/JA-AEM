clear; close all;
%% Parameters
addpath("Functions\");
addpath("Geometries\");
options.dmin = 0.3; %0.3e-3;
options.gamma = 1;
options.N = 100;
options.linear = 0;

%% Geometry test
tic
Geom = returnTwoSquares(options);
%%
[Geom] = linkNghbrBndrs(Geom);
Geom = returnNodes(Geom);
for nn = 1:length(Geom)
    [Geom(nn).G,Geom(nn).H] = calculateMatrixGHMultiBoundary(Geom(nn));
    [Geom(nn).G,Geom(nn).H,Geom(nn).bndr_val,Geom(nn).bndr_idx] = returnMultiDomainBndrVal(Geom(nn),Geom(nn).G,Geom(nn).H);
end
[A,B] = constructMultiDomain(Geom);
[Geom] = returnSolutionsMultiDomain(A,B,Geom);

figure
hold on 
for kk = 1:length(Geom)
    for nn = 1:length(Geom(kk).boundary)
        idx0 = Geom(kk).boundary(nn).bndr_idx == 0;
        idx1 = Geom(kk).boundary(nn).bndr_idx == 1;
        idx3 = Geom(kk).boundary(nn).bndr_idx == 3;

        plot(Geom(kk).boundary(nn).Xm(idx0),Geom(kk).boundary(nn).Ym(idx0),".k","MarkerSize",8)
        plot(Geom(kk).boundary(nn).Xm(idx1),Geom(kk).boundary(nn).Ym(idx1),".r","MarkerSize",8)
        plot(Geom(kk).boundary(nn).Xm(idx3),Geom(kk).boundary(nn).Ym(idx3),".m","MarkerSize",8)
        addArrowLine(Geom(kk).boundary(nn).xgeom,Geom(kk).boundary(nn).ygeom)    
    end
end
grid on 
xlabel("x-axis [a.u.]","Interpreter","latex","FontSize",12)
ylabel("y-axis [a.u.]","Interpreter","latex","FontSize",12)

%% Domain solutions
% Square 1
x1 = 0.001:0.03:4.999; 
y1 = 0.001:0.03:4.999;
[Xplot(1).val,Yplot(1).val] = meshgrid(x1,y1);

x2 = 5.001:0.03:9.999;
y2 = 0.001:0.03:4.999;
[Xplot(2).val,Yplot(2).val] = meshgrid(x2,y2);

x3 = 0.001:0.03:4.999;
y3 = 5.001:0.03:9.999;
[Xplot(3).val,Yplot(3).val] = meshgrid(x3,y3);

dplot = 10;
dumax = 100;

% maxdx = zeros(1,length(Geom));
% maxdy = zeros(1,length(Geom));
for nn = 1:length(Geom)
    [uplot(nn).val] = calculateDomainValuesMulti(Xplot(nn).val,Yplot(nn).val,Geom(nn));
    [uplot_dx(nn).val, uplot_dy(nn).val] = calculateDomainValuesDerivativeMulti(Xplot(nn).val,Yplot(nn).val,Geom(nn));
    uplot_dx(nn).val(abs(uplot_dx(nn).val)>dumax) = NaN;
    uplot_dy(nn).val(abs(uplot_dy(nn).val)>dumax) = NaN;
    % maxdx(nn) = max(uplot_dx(nn).val,[],"all");
    % maxdy(nn) = max(uplot_dy(nn).val,[],"all");
end

% maxval = max([maxdx,maxdy]);
% 
% for nn = 1:length(Geom)
%     uplot_dx(nn).val = uplot_dx(nn).val/maxval; 
%     uplot_dy(nn).val = uplot_dy(nn).val/maxval;
% end

toc
figure
hold on 

% contourf(XplotT,YplotT,UplotT,"FaceAlpha",0.5)
% quiver(XplotT(1:dplot:end,1:dplot:end),YplotT(1:dplot:end,1:dplot:end), ...
%     Uplot_dxT(1:dplot:end,1:dplot:end),Uplot_dyT(1:dplot:end,1:dplot:end),"LineWidth",1,"Color","k","AutoScale","on")

for nn = 1:length(Geom)
    contourf(Xplot(nn).val,Yplot(nn).val,uplot(nn).val,"FaceAlpha",0.5)
    quiver(Xplot(nn).val(1:dplot:end,1:dplot:end),Yplot(nn).val(1:dplot:end,1:dplot:end), ...
        uplot_dx(nn).val(1:dplot:end,1:dplot:end),uplot_dy(nn).val(1:dplot:end,1:dplot:end),"LineWidth",1,"Color","k","AutoScale","on")
end
colormap("summer")
c = colorbar;
c.Label.String = 'Magnitude [a.u]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";

grid on
xlabel("x-axis [a.u.]","Interpreter","latex","FontSize",12)
ylabel("y-axis [a.u.]","Interpreter","latex","FontSize",12)