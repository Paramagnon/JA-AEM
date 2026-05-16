function sigma = sigma_alu(T)
    %----- Data used for function fitting 
    % "Electrical resistivity of aluminium and manganese", 1984
    rho = [0.0001, 0.000102, 0.000109, 0.00014, 0.000192, 0.000345, 0.000748, ...
        0.00186, 0.00451, 0.0180, 0.047, 0.0955, 0.1618, 0.2439, 0.338, 0.440, ...
        1.003, 1.584, 2.155, 2.417, 2.650, 2.733, 3.875];        

    Tfit = [1,2,4,7,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,273,293, ...
        300,400]; % [K]
    
    % Resistivity in [Ohm.m]
    rho_out = interp1(Tfit,rho,T,"spline")*1e-8;
    sigma = 1./rho_out;
end