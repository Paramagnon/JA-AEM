function rho = rho_alu(T)   
    % Using fixed value for rho, design does not include shrinkage effects
    % So dimensions are how they would be during operation and small
    % assumption made here as calculted from dV = 3 alpha V0 dT worst-case
    % is around 1% 

    % Aluminum 6061-T6
    % Density [kg/m^3]
    rho = 2700*ones(size(T));     
end