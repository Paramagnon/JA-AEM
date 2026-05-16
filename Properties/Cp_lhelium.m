function Cp = Cp_lhelium(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=2.1768&THigh=5&TInc=0.1&Digits=5&ID=C7440597&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    Cp_fit = [29.531,25.267,23.194,22.352,22.230,22.545,23.141,23.929,24.857,...
        25.900,27.045,28.292,29.648,31.128,32.757,34.570,36.617,38.969,41.738,...
        45.099,49.354,51.798]; % [J/(g.K)]

    T_fit = [2.177,2.277,2.377,2.477,2.577,2.677,2.777,2.877,2.977,3.077,3.177,...
        3.277,3.377,3.477,3.577,3.677,3.777,3.877,3.977,4.077,4.177,4.2238]; % [K]

    % Heat capacity [J/(kg.K)]
    Cp = interp1(T_fit,Cp_fit,T,"spline")*1000; % [J/(kg.K)]

end