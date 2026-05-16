function f = currentPulsingCoil(Geom,par,t)
    Jref = par.Jref;
    N = par.NFour;
    k_n = @(n) n*pi/par.tau;
    
    %--- Motion profile
    theta = 2*pi*par.F*t;
    dx = 0;
    %---
    f = zeros(sum([Geom.Ldmn]),1);  
    x = Geom(1).domain.Xn;
    
    Nn = [-N:-1,1:N];       
    Jmov = zeros(length(theta),length(x));
    Jref = Jref*cos(theta);
    
    for tt = 1:length(theta)
        dxc = dx(tt); 
        for nn = 1:length(Nn)
            n = Nn(nn);
            k_nn = k_n(n);
        
            Jmovn = (-1i*Jref*exp(1i*pi*dxc/par.tau)/(2*n*pi)) * ( (exp(1i*k_nn*par.x2a) + exp(1i*k_nn*par.x3a) ...
                - exp(1i*k_nn*par.x1a) - exp(1i*k_nn*par.x4a)) + (exp(1i*k_nn*par.x2b) ...
                + exp(1i*k_nn*par.x3b) - exp(1i*k_nn*par.x1b) - exp(1i*k_nn*par.x4b))*exp(1i*2*pi/3) ...
                + (exp(1i*k_nn*par.x2c) + exp(1i*k_nn*par.x3c) - exp(1i*k_nn*par.x1c) ...
                - exp(1i*k_nn*par.x4c))*exp(1i*4*pi/3) );          

            Jmov(tt,:) = Jmov(tt,:) + Jmovn*exp(-1i*k_nn*(x-dxc));            
        end    
    end

    f(1:Geom(1).Ldmn) = real(Jmov');

end
