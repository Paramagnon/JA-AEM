clear; close all;
%% Parameters
p.Brem  = 1.5;                    % [T]
p.mu0 = 4*pi*1e-7; 


p.taup = 2e-3;
p.taum = 1e-3;
p.N = 500;

%% Magnetization


prefac = (4*Brem/(mu0*pi));






magnetization = 1;


%% Sine test
close all
theta = pi/3;
N = 500;

disp("Recurrence")
tic 
sinout = recurrenceSin(theta,N);
toc


disp("Sine direct")
tic 
sinout2 = zeros(1,N);
for nn = 1:N
    sinout2(nn) = sin(nn*theta);
end

toc

figure
plot(sinout-sinout2)

[mysn1,mysn2] = returnMysn(p.taum,p.taup,p);    
figure
plot(mysn1-mysn2)


%% Functions 
function [mysn1, msyn2] = returnMysn(taum,taup,p)    
    disp("reccurence")
    tic
    ndiv = 1./(1:p.N);
    fac1 = recurrenceSin(pi/2,p.N);
    fac2 = recurrenceSin(pi*taum/(2*taup),p.N);

    mysn1 = ndiv.*fac1.*fac2;
    toc

    disp("normal")
    tic
    msyn2 = zeros(1,p.N);    
    n = 1:p.N;
    mysn2 = (1./n).*sin(n*pi/2).*sin(n*pi*taum/(2*taup));
    toc
    % msyn = (1/n)*sin(n*pi/2)*sin(n*pi*taum/(2*taup))
end


function sinout = recurrenceSin(theta,N)
    sinout = zeros(1,N);
    sinout(1) = sin(1*theta);
    sinout(2) = sin(2*theta);
    fac1 = 2*cos(theta);

    for nn = 3:N
        sinout(nn) = fac1*sinout(nn-1) - sinout(nn-2);
    end
end











