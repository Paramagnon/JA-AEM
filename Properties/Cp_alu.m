function Cp = Cp_alu(T)
    % Implement temperature dependent thermal conductivity
    % https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061

    a = 46.6467;
    b = -314.292;
    c = 866.662;
    d = -1298.3;
    e = 1162.27;
    f = -637.795;
    g = 210.351;
    h = -38.3094;
    k = 2.96344;
     
    % Thermal conductivity [W/(m.K)]
    T10 = log10(T);
    Cp = 10.^(a+b.*T10 + c.*T10.^2 + d.*T10.^3 + e.*T10.^4 + f.*T10.^5 + g.*T10.^6 + h.*T10.^7 + k.*T10.^8);
end