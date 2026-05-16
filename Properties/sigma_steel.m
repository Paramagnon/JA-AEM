function sigma = sigma_steel(T)
    %----- Data used for function fitting 
    % "Update of thermal conductivity and electrical resistivity of electrolytic 
    % iron, tungsten, and stainless steel", 1984
    rho = [593, 593, 593, 594, 594, 594, 594, 594, 594, 594, 594, 593, 593, ... 
        593, 592, 592, 593, 595, 597, 599, 606, 613, 622, 630, 639, 683, 724, ...
        767, 810, 885]; % [n.Ohm.m]

    Tfit = [2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,...
        100,150,200,250,300,400]; % [K]
    
    % Resistivity in [Ohm.m]
    rho_out = interp1(Tfit,rho,T,"spline")*1e-9;
    sigma = 1./rho_out;

end