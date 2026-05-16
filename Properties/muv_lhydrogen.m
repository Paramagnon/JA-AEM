function muv = muv_lhydrogen(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=13.957&THigh=25&TInc=0.25&Digits=5&ID=C1333740&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    muv_fit = [25.758,24.956,24.199,23.484,22.806,22.163,21.551,20.967,20.411,...
        19.879,19.370,18.882,18.414,17.965,17.533,17.117,16.717,16.331,15.959,...
        15.600,15.253,14.917,14.592,14.278,13.973,13.677,13.490]; % [uPa.s]

    T_fit = [13.957,14.207,14.457,14.707,14.957,15.207,15.457,15.707,15.957,...
        16.207,16.457,16.707,16.957,17.207,17.457,17.707,17.957,18.207,18.457,...
        18.707,18.957,19.207,19.457,19.707,19.957,20.207,20.369]; % [K]

    % Dynamic viscosity [Pa.s]
    muv = interp1(T_fit,muv_fit,T,"spline")*1e-6; % [Pa.s]
   
end