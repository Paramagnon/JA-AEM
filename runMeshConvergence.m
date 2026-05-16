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

options.SC = 0;
options.Ec = 1e-4;
options.n = 30;
options.Jc = 1.149e10;
options.mu_sc = 4*pi*1e-7;
options.d_tape = 4e-6;
options.gamma = 10; 
options.Iset = 10;
options.trise = 0.1; 
options.ACstart = 0.02;
options.ACrise = 0.001;

options.maxit = 500;
options.maxerror = 1e-6;

options.file = 'test.gif';
options.debug = 1;                  % Display time step number
options.runGIF = 0;                 % Run GIF output of transient
options.dlplot = 2e6;               % Plot density for contourf
options.crcl = 0;                   % Plot black circles for the water pipes
options.pltL = 0;                   % Plot only internal nodes
options.pltlines = 1;

options.dt = 0.01;
options.tmax = 0.1;
%% Mesh settings
[par] = loadParCoilSolver();
par.T0 = 169;                                % [K]
options.dmin = 1e-4;
options.gamma = 5;
options.linear = 1;
options.Nmin = 5;

%% Setup - Coil solver
N = 600:100:2000;
dL = 0.1e6:0.1e6:2.5e6;

Linf2 = zeros(length(N),length(dL)); 
L22 = zeros(length(N),length(dL));

par.Jref = 10e6;
par.F = 5;
options.par = par;
options.Iset = 300;

for nn = 1:length(N)
    options.N = N(nn);
    for ll = 1:length(dL)
        options.dL = dL(ll);
        
        options.dLx = sqrt(options.dL);
        options.dLy = sqrt(options.dL);
        options.Nx = options.N;
        options.Ny = options.N;
        options.dy = par.h_turn/(par.h_turn*2*options.Ny-1);

        [Geom] = returnCoilSolver(options,par);
        [Geom] = linkNghbrBndrs(Geom);
        [Geom] = returnNodesAEM(Geom);
        [SimGeom] = runSimInit(Geom,par,options);

        [SimGeom] = runTransient(Geom, SimGeom, options, @currentWaveform);

        Linf2(nn,ll) = max(SimGeom.SimSol(end).uL);
        L22(nn,ll) = sum(sqrt(SimGeom.SimSol(end).uL.^2))/length(SimGeom.SimSol(end).uL);
        disp(['Case: ' num2str((nn-1)*length(dL)+ll) ', L2: ' num2str(L22(nn,ll)) ', Linf: ' num2str(Linf2(nn,ll))])
    end
end

%% Figures
close all
[DL,DN] = meshgrid(dL(2:end),N(3:end));
% DN = repmat(N(2:end),length(dL),1);
% DL = repmat(dL,1,length(N)-1);
% DN = reshape(DN,[],1);

Z2 = L2(3:end,2:end)/L2(end,end);
Zinf = Linf(3:end,2:end)/Linf(end,end);


% Zcolor = returnZColor(Zinf);

figure
contourf(DN,DL,Z2)
c = colorbar;
c.Label.String = 'Normalized $L_2$ [-]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";
xlabel("Line density N [-]","fontsize",12,"interpreter","latex")
ylabel("Domain density dL [-]","fontsize",12,"interpreter","latex")
colorMapAltair

figure
contourf(DN,DL,Zinf)
c = colorbar;
c.Label.String = 'Normalized $L_{\infty}$ [-]';
c.Label.FontSize = 12;
c.Label.Interpreter = "latex";
xlabel("Line density N [-]","fontsize",12,"interpreter","latex")
ylabel("Domain density dL [-]","fontsize",12,"interpreter","latex")
colorMapAltair

%% Functions

function Zcolor = returnZColor(Z)
    N = size(Z,2);
    Zcolor = zeros(size(Z,1)*N,3);
    
    colorMapAltairRGB = [1, 3, 114; ...
                      4, 1, 161; ...
                      5, 2, 192; ...
                      79, 5, 235; ...
                      135, 2, 206; ...
                      175, 3, 177; ...
                      206, 5, 139; ...
                      235, 4, 77; ...
                      251, 1, 81; ...
                      253, 72, 10; ...
                      253, 122, 7; ...
                      250, 162, 3; ...
                      252, 192, 3; ...
                      251, 221, 1; ...
                      255, 237, 9; ...
                      252, 253, 7; ...
                      251, 254, 178; ...
                     ]/255;

    Zmin = min(Z,[],"all");
    Zmax = max(Z,[],"all");
    rgbidx = linspace(Zmin,Zmax,17);

    for ii = 1:size(Z,1)
        
        r = interp1(rgbidx,colorMapAltairRGB(:,1),Z(ii,:));
        g = interp1(rgbidx,colorMapAltairRGB(:,2),Z(ii,:));
        b = interp1(rgbidx,colorMapAltairRGB(:,3),Z(ii,:));
    
        Zcolor((ii-1)*N+1:ii*N,:) = [r' g' b'];

    end
end