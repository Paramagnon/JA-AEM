function kappa = kappa_lnitrogen(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=63.151&THigh=80&TInc=0.5&Digits=5&ID=C7727379&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    kappa_fit = [0.17335,0.17234,0.17133,0.17033,0.16932,0.16831,0.16731,0.16630,...
        0.16530,0.16429,0.16329,0.16229,0.16128,0.16028,0.15928,0.15827,0.15727,...
        0.15627,0.15527,0.15427,0.15327,0.15226,0.15126,0.15026,0.14926,0.14826,...
        0.14726,0.14626,0.14526,0.14485]; % [W/(m.K)]

    T_fit = [63.151,63.651,64.151,64.651,65.151,65.651,66.151,66.651,67.151,...
        67.651,68.151,68.651,69.151,69.651,70.151,70.651,71.151,71.651,72.151,...
        72.651,73.151,73.651,74.151,74.651,75.151,75.651,76.151,76.651,77.151,77.355]; % [K]

    % Thermal conductivity [W/(m.K)]   
    kappa = interp1(T_fit,kappa_fit,T,"spline"); % [W/(m.K)]

end