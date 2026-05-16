function sigma = sigma_copper(T)
    %----- Data used for function fitting 
    % "Electrical resistivity of copper, gold, palladium, and silver", 1979
    rho = [0.002, 0.002, 0.002, 0.00202, 0.00218, 0.0028, 0.00449, 0.00828, ...
        0.0147, 0.0239, 0.0358, 0.0518, 0.0727, 0.0971, 0.154, 0.215, 0.281, ...
        0.348, 0.522, 0.699, 0.874, 1.046, 1.217, 1.387, 1.543, 1.678, 1.725, ...
        2.063]; 

    Tfit = [1,4,7,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,...
        225,250,273.15,293,300,350]; % [K]
    
    % Resistivity in [Ohm.m]
    rho_out = interp1(Tfit,rho,T,"spline")*1e-8;
    sigma = 1./rho_out;

end