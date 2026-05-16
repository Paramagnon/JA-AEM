par.NFour = 300;
par.tau_m = 30e-3;
par.tau = 60e-3;
tau_c = 40e-3;

dc = 5e-3;
wc = 0.012;
id = 0.011;

offs = wc;
par.x1a = dc/2+offs;
par.x2a = dc/2+wc+offs;
par.x3a = dc/2+wc+id+offs;
par.x4a = dc/2+wc+id+wc+offs;

par.x1b = dc/2+tau_c+offs;
par.x2b = dc/2+wc+tau_c+offs;
par.x3b = dc/2+wc+id+tau_c+offs;
par.x4b = dc/2+wc+id+wc+tau_c+offs;

par.x1c = dc/2+2*tau_c+offs;
par.x2c = dc/2+wc+2*tau_c+offs;
par.x3c = dc/2+wc+id+2*tau_c+offs;
par.x4c = dc/2+wc+id+wc+2*tau_c+offs;

par.F = 20;
options.par = par;

[an, bn] = boundaryFourierInnit(par);
xplt = linspace(0,0.12,200);



t = 0.032;
options.t = t;
[jplt] = boundaryFourierCurrent(an,bn,xplt,options,t);
figure(1)
plot(xplt,jplt)

% for ii = 1:30
%     t = t+0.0001;
%     options.t = t;
%     [jplt] = boundaryFourierCurrent(an,bn,xplt,options,t);
%     figure(1)
%     plot(xplt,jplt)
%     pause(1)
% end



%% Functions

function [jplt] = boundaryFourierCurrent(an,bn,xplt,options,t)
    par = options.par;
    options.t = t;
    Jref = currentAC(options, par);    
    k_n = @(n) n*pi/par.tau;
    
    %--- Motion profile
    theta = 2*pi*par.F*(t-options.ACstart);
    tau_m = 0.03;
    x = -4*tau_m*sin(theta);
    theta = pi*x/tau_m;
    %--- Coil currents
    Ja = cos(theta+pi)*Jref
    Jb = cos(theta+pi-2*pi/3)*Jref %--> equal Jc in comsol
    Jc = cos(theta+pi+2*pi/3)*Jref %--> equal Jb in comsol

    %--- Actual values
    Nn = [1:par.NFour]';  
    k_nn = k_n(Nn);    

    % theta = 0; %---> Turn off movement, if you want pulsing coils
    costheta = cos(k_nn'*theta);
    sintheta = sin(k_nn'*theta);
    
    an_nn = Ja*(an(1,:).*costheta + bn(1,:).*sintheta) + Jb*(an(2,:).*costheta + bn(2,:).*sintheta) + Jc*(an(3,:).*costheta + bn(3,:).*sintheta);
    bn_nn = Ja*(bn(1,:).*costheta - an(1,:).*sintheta) + Jb*(bn(2,:).*costheta - an(2,:).*sintheta) + Jc*(bn(3,:).*costheta - an(3,:).*sintheta);

    jplt = an_nn*sin(k_nn*xplt) + bn_nn*cos(k_nn*xplt);
end

function Jref = currentAC(options, par)  
    if options.t < options.ACstart
        Jref = 0;
    elseif options.t < options.ACstart + options.ACrise
        Jref = (options.t - options.ACstart)*par.Jref/options.ACrise;
    else
        Jref = par.Jref; 
    end
end

function [an, bn] = boundaryFourierInnit(par)  
    k_n = @(n) n*pi/par.tau; 
    Nn = [1:par.NFour];   
    
    an = zeros(3,length(Nn)); % Phase [a,b,c] 
    bn = zeros(3,length(Nn)); % Phase [a,b,c]
     
    k_nn = k_n(Nn);

    an(1,:) = (1./(k_nn*par.tau)).*(-cos(k_nn*par.x2a) + cos(k_nn*par.x1a) + cos(k_nn*par.x4a) - cos(k_nn*par.x3a));
    an(2,:) = (1./(k_nn*par.tau)).*(-cos(k_nn*par.x2b) + cos(k_nn*par.x1b) + cos(k_nn*par.x4b) - cos(k_nn*par.x3b));
    an(3,:) = (1./(k_nn*par.tau)).*(-cos(k_nn*par.x2c) + cos(k_nn*par.x1c) + cos(k_nn*par.x4c) - cos(k_nn*par.x3c));

    bn(1,:) = (1./(k_nn*par.tau)).*(sin(k_nn*par.x2a) - sin(k_nn*par.x1a) - sin(k_nn*par.x4a) + sin(k_nn*par.x3a));
    bn(2,:) = (1./(k_nn*par.tau)).*(sin(k_nn*par.x2b) - sin(k_nn*par.x1b) - sin(k_nn*par.x4b) + sin(k_nn*par.x3b));
    bn(3,:) = (1./(k_nn*par.tau)).*(sin(k_nn*par.x2c) - sin(k_nn*par.x1c) - sin(k_nn*par.x4c) + sin(k_nn*par.x3c));
end