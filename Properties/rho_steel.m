function rho = rho_steel(T)   
    % Using fixed value for rho, design does not include shrinkage effects
    % So dimensions are how they would be during operation and small
    % assumption made here as calculted from dV = 3 alpha V0 dT worst-case
    % is around 1% 

    % 304 Stainless
    % Density [kg/m^3]
    rho = 7930*ones(size(T));     
end