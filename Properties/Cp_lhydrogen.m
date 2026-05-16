function Cp = Cp_lhydrogen(T)  
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=13.957&THigh=25&TInc=0.25&Digits=5&ID=C1333740&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    Cp_fit = [70.101,70.720,71.440,72.236,73.091,73.991,74.926,75.889,76.876,...
        77.883,78.909,79.951,81.011,82.088,83.182,84.296,85.429,86.584,87.762,...
        88.965,90.193,91.450,92.738,94.057,95.411,96.802,97.725]; % [J/(g.K)]

    T_fit = [13.957,14.207,14.457,14.707,14.957,15.207,15.457,15.707,15.957,...
        16.207,16.457,16.707,16.957,17.207,17.457,17.707,17.957,18.207,18.457,...
        18.707,18.957,19.207,19.457,19.707,19.957,20.207,20.369]; % [K]

    % Specific heat [J/(kg.K)]
    Cp = interp1(T_fit,Cp_fit,T,"spline")*1000; % [J/(kg.K)]
       
end