function kappa = kappa_lhelium(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=2.1768&THigh=5&TInc=0.1&Digits=5&ID=C7440597&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    kappa_fit = [0.013718,0.014191,0.014629,0.015037,0.015420,0.015781,0.016122,...
        0.016443,0.016745,0.017028,0.017292,0.017536,0.017758,0.017959,0.018138,...
        0.018292,0.018420,0.018520,0.018591,0.018629,0.018631,0.018619]; % [W/(m.K)]

    T_fit = [2.177,2.277,2.377,2.477,2.577,2.677,2.777,2.877,2.977,3.077,3.177,...
        3.277,3.377,3.477,3.577,3.677,3.777,3.877,3.977,4.077,4.177,4.2238]; % [K]

    % Thermal conductivity [W/(m.K)]   
    kappa = interp1(T_fit,kappa_fit,T,"spline"); % [W/(m.K)]

end