function Cp = Cp_lnitrogen(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=63.151&THigh=80&TInc=0.5&Digits=5&ID=C7727379&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    Cp_fit = [19.997,20.006,20.014,20.023,20.032,20.041,20.051,20.061,20.071,...
        20.082,20.093,20.105,20.118,20.131,20.144,20.158,20.173,20.188,20.204,...
        20.221,20.238,20.256,20.275,20.295,20.315,20.336,20.358,20.381,20.405,20.415]; % [J/(g.K)]

    T_fit = [63.151,63.651,64.151,64.651,65.151,65.651,66.151,66.651,67.151,...
        67.651,68.151,68.651,69.151,69.651,70.151,70.651,71.151,71.651,72.151,...
        72.651,73.151,73.651,74.151,74.651,75.151,75.651,76.151,76.651,77.151,77.355]; % [K]

    % Heat capacity [J/(kg.K)]
    Cp = interp1(T_fit,Cp_fit,T,"spline")*1000; % [J/(kg.K)]
   
end