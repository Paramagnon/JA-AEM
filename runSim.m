clear; close all;
%% Parameters
addpath("Functions");
addpath("AEM\Matrices\");
addpath("AEM\Mesh\");
addpath("AEM\TimeSolver\");
addpath("AEM\Properties\");
addpath("AEM\Parameters\");
addpath("AEM\Plotting\");
addpath("AEM\Geometries\");

options.dmin = 1e-4;
options.gamma = 5;
options.N = 1000;
options.Nx = options.N;
options.Ny = options.N;
options.Nmin = 3;
options.linear = 1;
options.dL = 0.6e6;
options.dLx = sqrt(options.dL);
options.dLy = sqrt(options.dL);
options.dlplot = 2e6;
options.pltlines = 1;

options.SC = 1;                     % Toggle SC properties on/off
options.Ec = 1e-4;
options.n = 30;
options.Jc = 1.149e10;
options.mu_sc = 4*pi*1e-7;
options.d_tape = 4e-6;
options.gamma = 10; 
options.Iset = 10;
options.trise = 0.1; 

options.maxit = 500;
options.maxerror = 1e-6;
%% Setup - Coil solver
tic
[par] = loadParCoilSolver();
par.T0 = 169;                                % [K]
par.NFour = 100;
par.srcgeom = 1;
[Geom] = returnCoilSolver(options,par);
% [Geom] = returnContinuousTest(options);
[Geom] = linkNghbrBndrs(Geom);
[Geom] = returnNodesAEM(Geom);
[SimGeom] = runSimInit(Geom,par,options);
toc
%% Show Mesh
plotMesh(SimGeom.PlotGeom,0)
plotMesh(Geom,0)
%% Steady state solution
[SimGeom] = runSteadyState(Geom, SimGeom, options);
options.pltL = 1;

figure
[~] = plotContourf(Geom, SimGeom.PlotGeom, SimGeom.SimSol, options, 1);
c = colorbar;
c.Label.String = 'Magnitude [Wb/m]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";

grid on
xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
axis square
%% Transient solution
tic
options.file = 'SC_full2.gif';
options.dt = 0.005;
options.tmax = 0.08;
options.debug = 1;                  % Display time step number
options.runGIF = 0;                 % Run GIF output of transient
options.dlplot = 2e6;               % Plot density for contourf
options.pltL = 0;                   % Plot only internal nodes
options.dy = par.h_turn/(par.h_turn*2*options.Ny-1);

par.Jref = 0e6;
par.F = 5;
options.par = par;
options.Iset = 300;
options.trise = 0.01; 
options.ACstart = 0.02;
options.ACrise = 0.001;

[SimGeom] = runTransient(Geom, SimGeom, options, @currentWaveform);
% [SimGeom] = runTransient(Geom, SimGeom, options, @currentSine);
toc
%% Plot output
plotGIF(Geom, SimGeom, options);   
% plotSCCurrents(SimGeom,options);