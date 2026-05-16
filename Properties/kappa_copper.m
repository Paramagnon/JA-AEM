function kappa = kappa_copper(T,RRR)
    % Implement temperature dependent thermal conductivity
    % https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
    
    switch RRR
        case 50 
            %--- RRR = 50
            % Range: 4-300K
            a = 1.8743;
            b = -0.41538;
            c = -0.6018;
            d = 0.13294; 
            e = 0.26426; 
            f = -0.0219;
            g = -0.051276;
            h = 0.0014871;
            k = 0.003723;
        case 300 
            %--- RRR = 300
            % Range: 4-300K
            a = 1.357;
            b = 0.3981;
            c = 2.669;
            d = -0.1346; 
            e = -0.6683; 
            f = 0.01342;
            g = 0.05773;
            h = 0.0002147;
            k = 0;
        otherwise
            disp('No matching copper material')
    end
    
    kappa = 10.^(( a + c*T.^0.5 + e*T + g*T.^1.5 + k*T.^2 ) ./ ( 1+ b*sqrt(T) + d*T + f*T.^1.5 + h*T.^2));
end