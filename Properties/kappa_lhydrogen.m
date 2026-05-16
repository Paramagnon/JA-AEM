function kappa = kappa_lhydrogen(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=13.957&THigh=25&TInc=0.25&Digits=5&ID=C1333740&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    kappa_fit = [0.097630,0.098046,0.098448,0.098836,0.099210,0.099569,0.099915,...
        0.10025,0.10056,0.10087,0.10115,0.10143,0.10168,0.10193,0.10216,0.10237,...
        0.10257,0.10275,0.10292,0.10308,0.10322,0.10334,0.10345,0.10354,0.10362,...
        0.10368,0.10371]; % [W/(m.K)]

    T_fit = [13.957,14.207,14.457,14.707,14.957,15.207,15.457,15.707,15.957,...
        16.207,16.457,16.707,16.957,17.207,17.457,17.707,17.957,18.207,18.457,...
        18.707,18.957,19.207,19.457,19.707,19.957,20.207,20.369]; % [K]

    % Thermal conductivity [W/(m.K)]   
    kappa = interp1(T_fit,kappa_fit,T,"spline"); % [W/(m.K)]
    
end