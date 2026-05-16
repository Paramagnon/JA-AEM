% clear; close all;
%% Add paths
ogfolder = pwd;
cd('..\');
addpath("Parameters\");
addpath("TimeSolver\");
cd(ogfolder);
%% Load variables
% [par] = loadParCoilSolver();
par = loadParCoilCOMSOL();
Geom(1).domain.Xn = linspace(par.x0,par.n_m*par.tau_m,500);
Geom(1).Ldmn = 500;
par.Jref = 10e6;
par.F = 100;
par.NFour = 1000;
t = 0:0.0005:0.080;

figure("Position",[1,1,735,164],"Units","pixels")
% figure
p = plot(nan,nan,"LineWidth",2);
p.XData = Geom(1).domain.Xn;


grid on
ylim([-1.2*par.Jref, 1.2*par.Jref])
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("Current [A/m$^2$]","Interpreter","latex","FontSize",12)
% axis square 

disp(['Timestep: ' num2str(1)])
p.YData = currentPulsingCoil(Geom,par,t(1));
exportgraphics(gcf,'CurrentWaveform.gif');

for tt = 2:length(t)   
    disp(['Timestep: ' num2str(tt)])
    p.YData = currentPulsingCoil(Geom,par,t(tt));
    exportgraphics(gcf,'CurrentWaveform.gif','Append',true);
end