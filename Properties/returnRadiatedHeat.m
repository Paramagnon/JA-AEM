function qn = returnRadiatedHeat(Top,emiss)
    % Linearized radiated heat, for a difference in temperature of 220 K
    % the error for the linearized radiated heat is 0.3% 
    SB = 5.67e-8; % Boltzmann constant    
    qn = SB*emiss*Top.^4;
end