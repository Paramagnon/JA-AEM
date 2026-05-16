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