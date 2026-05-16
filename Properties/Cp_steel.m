function Cp = Cp_steel(T)
    % Implement temperature dependent thermal conductivity
    % https://trc.nist.gov/cryogenics/materials/304Stainless/304Stainless_rev.htm 

    %----- 304 Stainless (UNS S30400) -----
    % https://trc.nist.gov/cryogenics/materials/310%20Stainless/310Stainless_rev.htm
    % Range: 4-47K
    a = 22.0061;
    b = -127.5528;
    c = 303.647;
    d = -381.0098;
    e = 274.0328;
    f = -112.9212;
    g = 24.7593;
    h = -2.239153;
    k = 0;
    
    % Heat capacity [J/(kg.K)]
    T10 = log10(T);
    Cp = 10.^(a+b*T10 + c*T10.^2 + d*T10.^3 + e*T10.^4 + f*T10.^5 + g*T10.^6 + h*T10.^7 + k*T10.^8);
end