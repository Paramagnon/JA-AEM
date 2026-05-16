function rho = rho_lhelium(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=2.1768&THigh=5&TInc=0.1&Digits=5&ID=C7440597&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    rho_fit = [147.92,147.71,147.35,146.88,146.32,145.67,144.96,144.17,143.32,...
        142.39,141.40,140.34,139.20,137.97,136.66,135.23,133.70,132.02,130.18,...
        128.15,125.85,124.67]; % [kg/m^3]

    T_fit = [2.177,2.277,2.377,2.477,2.577,2.677,2.777,2.877,2.977,3.077,3.177,...
        3.277,3.377,3.477,3.577,3.677,3.777,3.877,3.977,4.077,4.177,4.2238]; % [K]

    % Density [kg/m^3]   
    rho = interp1(T_fit,rho_fit,T,"spline"); % [kg/m^3]
    
end