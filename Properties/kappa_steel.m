function kappa = kappa_steel(T,ST)
    % Implement temperature dependent thermal conductivity
    % https://trc.nist.gov/cryogenics/materials/materialproperties.htm
    
    %----- 310 Stainless (UNS S31000) -----
    % https://trc.nist.gov/cryogenics/materials/310%20Stainless/310Stainless_rev.htm
    % Range: 4-300K
    a = -0.81907;
    b = -2.1967;
    c = 9.1059;
    d = -13.078;
    e = 10.853;
    f = -5.1269;
    g = 1.2583;
    h = -0.12395;
    k = 0;

    % Thermal conductivity [W/(m.K)]
    T10 = log10(T);
    kappa = 10.^(a+b*T10 + c*T10.^2 + d*T10.^3 + e*T10.^4 + f*T10.^5 + g*T10.^6 + h*T10.^7 + k*T10.^8);
end