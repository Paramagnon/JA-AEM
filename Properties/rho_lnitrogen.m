function rho = rho_lnitrogen(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=63.151&THigh=80&TInc=0.5&Digits=5&ID=C7727379&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    rho_fit = [867.37,865.32,863.26,861.19,859.12,857.04,854.95,852.86,850.76,...
        848.65,846.54,844.41,842.28,840.14,838.00,835.84,833.68,831.51,829.33,...
        827.14,824.94,822.73,820.52,818.29,816.06,813.81,811.56,809.29,807.02,806.08]; % [kg/m^3]

    T_fit = [63.151,63.651,64.151,64.651,65.151,65.651,66.151,66.651,67.151,...
        67.651,68.151,68.651,69.151,69.651,70.151,70.651,71.151,71.651,72.151,...
        72.651,73.151,73.651,74.151,74.651,75.151,75.651,76.151,76.651,77.151,77.355]; % [K]

    % Density [kg/m^3]   
    rho = interp1(T_fit,rho_fit,T,"spline"); % [kg/m^3]
   
end