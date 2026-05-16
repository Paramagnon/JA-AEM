function rho = rho_lhydrogen(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=13.957&THigh=25&TInc=0.25&Digits=5&ID=C1333740&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    rho_fit = [77.083,76.879,76.672,76.464,76.253,76.039,75.822,75.602,75.380,...
        75.154,74.925,74.693,74.458,74.218,73.976,73.729,73.479,73.225,72.967,...
        72.705,72.439,72.168,71.893,71.613,71.328,71.039,70.848]; % [kg/m^3]

    T_fit = [13.957,14.207,14.457,14.707,14.957,15.207,15.457,15.707,15.957,...
        16.207,16.457,16.707,16.957,17.207,17.457,17.707,17.957,18.207,18.457,...
        18.707,18.957,19.207,19.457,19.707,19.957,20.207,20.369]; % [K]

    % Density [kg/m^3]   
    rho = interp1(T_fit,rho_fit,T,"spline"); % [kg/m^3]
    
end