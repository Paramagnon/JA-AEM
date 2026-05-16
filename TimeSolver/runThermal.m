function [alpha,SimGeom] = runThermal(SimGeom,alpha_prev,options) 
    N = size(SimGeom.H00T,2);
    %---------------- Matrices ---------------------
    dtU = options.dt*SimGeom.UxxUyyT;
    if options.Euler
        % Backward Euler:
        dtf = options.dt*SimGeom.fT;
    else
        % Crank–Nicolson:
        dtf = options.dt*(SimGeom.fT+SimGeom.fT_prev);
    end
    C0 = (SimGeom.rho.*SimGeom.Cp./SimGeom.kappa);
    CS = C0.*SimGeom.S00T;
    dC3 = zeros(2*N,1);
    dC3(N+1:2*N) = SimGeom.C3T - SimGeom.C3T_prev;
    dz = ([SimGeom.H00T, -SimGeom.G00T]*SimGeom.invKT)*dC3;

    %------------ Next time step --------- 
    if options.Tupdate  
        if options.Euler
            % Backward Euler:
            alpha = (dtU-CS)\(dtf + C0.*dz - CS*alpha_prev);
        else            
            % Crank–Nicolson:
            alpha = (dtU-2*CS)\(dtf + 2*C0.*dz - (2*CS + dtU)*alpha_prev);
        end 
    else
        if options.Euler
            % Backward Euler:
            alpha = SimGeom.invdtUCST*(dtf + C0.*dz - CS*alpha_prev);       
        else
            % Crank–Nicolson:
            alpha = SimGeom.invdtUCST*(dtf + 2*C0.*dz - (2*CS + dtU)*alpha_prev);
        end
    end
    [SimGeom.eT,SimGeom.z00T] = updateC3T(SimGeom);
    
    SimGeom.TLT = SimGeom.S00T*alpha + SimGeom.z00T;
end