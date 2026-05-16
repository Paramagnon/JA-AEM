function muv = muv_lnitrogen(T)
    %----- Data used for function fitting  (1 Atm)
    % https://webbook.nist.gov/cgi/fluid.cgi?P=1&TLow=63.151&THigh=80&TInc=0.5&Digits=5&ID=C7727379&Action=Load&Type=IsoBar&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

    muv_fit = [311.95,303.52,295.43,287.64,280.15,272.95,266.02,259.34,252.91,...
        246.71,240.74,234.98,229.42,224.06,218.88,213.88,209.05,204.39,199.88,...
        195.51,191.29,187.21,183.25,179.42,175.71,172.12,168.63,165.25,161.97,160.66]; % [uPa.s]

    T_fit = [63.151,63.651,64.151,64.651,65.151,65.651,66.151,66.651,67.151,...
        67.651,68.151,68.651,69.151,69.651,70.151,70.651,71.151,71.651,72.151,...
        72.651,73.151,73.651,74.151,74.651,75.151,75.651,76.151,76.651,77.151,77.355]; % [K]

    % Dynamic viscosity [Pa.s]
    muv = interp1(T_fit,muv_fit,T,"spline")*1e-6; % [Pa.s]

end