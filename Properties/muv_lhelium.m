function muv = muv_lhelium(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=2.1768&THigh=5&TInc=0.1&Digits=5&ID=C7440597&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    muv_fit = [37.766,38.530,38.990,39.222,39.282,39.212,39.041,38.790,38.477,...
        38.112,37.705,37.263,36.789,36.286,35.757,35.202,34.619,34.006,33.359,...
        32.670,31.927,31.555]; % [uPa.s]

    T_fit = [2.177,2.277,2.377,2.477,2.577,2.677,2.777,2.877,2.977,3.077,3.177,...
        3.277,3.377,3.477,3.577,3.677,3.777,3.877,3.977,4.077,4.177,4.2238]; % [K]

    % Dynamic viscosity [Pa.s]
    muv = interp1(T_fit,muv_fit,T,"spline")*1e-6; % [Pa.s]

end