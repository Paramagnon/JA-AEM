clear; close all;
%% Paths
addpath("Functions");
addpath("AEM\Matrices\");
addpath("AEM\Mesh\");
addpath("AEM\TimeSolver\");
addpath("AEM\Properties\");
addpath("AEM\Parameters\");
addpath("AEM\Plotting\");
addpath("AEM\Geometries\");
%% Parameters
options.dmin = 1e-4; %1e-4;
options.gamma = 0.6; % Smaller is more samples
options.N = 2500; %1000;
options.Nx = 1500; %500;
options.Ny = options.N;
options.Nmin = 6;
options.linear = 1;
options.dL = 0.5e6; %1.5e6;
options.dLx = 750;
options.dLy = 3200;
options.dlplot = 1.5e6;
options.pltlines = 1;

options.SC = 1;                 % Toggle SC properties on/off
options.Ec = 1e-4;
options.n = 30;
options.Jc = 1.149e10;
options.mu_sc = 4*pi*1e-7;
options.d_tape = 2e-6;
options.Iset = 10;
% options.Ifrac = 0.8;          % Fraction of Ic
options.trise = 0.01; 
options.ACstart = 0.02;
options.ACrise = 0.001;

options.maxit = 20;
options.maxerror = 1e-6;
options.Tupdate = 0;            % Update properties every temp step
options.Euler = 1;              % Use Backward Euler (1) or Crank-Nicholson (0)

%% Setup Geometry
tic
par = loadParCoilCOMSOL();
par.T0 = 293;                               % [K]
% options.Jc = returnJcHTS(par.T0,2);         % (T0 [K], B0 [T]) --> Returns zero if above critical temperature
% options.Iset = options.Jc*par.h_turn*options.d_tape*options.Ifrac;
par.srcgeom = 2;
options.NySC = (30+1)/par.h_turn;
[Geom] = returnCoilCOMSOL(options,par);
[SimGeom] = runSimInit(Geom,par,options);
toc
%% Show Mesh
% plotMesh(SimGeom.PlotGeom,0)
plotMesh(Geom,0)
plotMeshThermal(Geom,0)

%% Transient solution
tic
options.file = 'COMSOL.gif';
options.dt = 1e-4;
options.tmax = 0.03;
options.debug = 1;                  % Display time step number
options.runGIF = 0;                 % Run GIF output of transient
options.pltL = 0;                   % Plot only internal nodes
% options.dy = par.h_turn/(par.h_turn*options.NySC-1);

par.Jref = 10e6;
par.F = 100;
options.par = par;
options.Iset = 73.54;
[SimGeom] = runTransient(Geom, SimGeom, options, @currentWaveform);
SimGeom.runidx = 100;
toc

%% 
% plotGIF(Geom, SimGeom, options);   
% plotGIFThermal(Geom, SimGeom, options);   
plotSCCurrents(SimGeom,options, 0);