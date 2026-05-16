function kappa = kappa_alu(T,ALU)
    % Implement temperature dependent thermal conductivity
    % https://trc.nist.gov/cryogenics/materials/materialproperties.htm
 
    switch ALU
        case 1100 
            %----- 1100 Aluminum (UNS A91100) -----
            % https://trc.nist.gov/cryogenics/materials/1100%20Aluminum/1100%20Aluminum_rev.htm
            % Range: 4-300K
            a = 23.39172;
            b = -148.5733;
            c = 422.1917;
            d = -653.6664;
            e = 607.0402;
            f = -346.152;
            g = 118.4276;
            h = -22.2781;
            k = 1.770187;
        case 3003
            %----- 3003F Aluminum (UNS A93003) -----
            % https://trc.nist.gov/cryogenics/materials/3003F%20Aluminum/3003FAluminum_rev.htm
            % Range: 4-300K    
            a = 0.63736;
            b = -1.1437;
            c = 7.4624;
            d = -12.6905;
            e = 11.9165;
            f = -6.18721;
            g = 1.63939;
            h = -0.172667;
            k = 0;
        otherwise   
            disp('No matching aluminium material')
    end
    % Thermal conductivity [W/(m.K)]
    T10 = log10(T);
    kappa = 10.^(a+b*T10 + c*T10.^2 + d*T10.^3 + e*T10.^4 + f*T10.^5 + g*T10.^6 + h*T10.^7 + k*T10.^8);
end