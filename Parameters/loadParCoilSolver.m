function par = loadParCoilSolver()
    par.x0 = 0;
    par.y0 = 0;
        
    par.h_c = 10e-3;                    % Coil height/2 (symmetry in z-axis) [m]
    par.h_g = 3e-3;                      % Air gap height [m]
    par.h_cw = 5e-3;                    % Cryo wall height [m]
    par.h_turn = 8e-3;                      % Height per turn [m]
    par.h_p = 5e-3;                
    
    par.tau_m = 15e-3;                  % Magnet period [m] 
    par.w_turn = 114e-6;                    % Width per turn [m]
    
    par.r_pipes = 1e-3;                 % Radius of pipes [m]

    par.n_m = 2;                            % Number of magnets [-]
    par.n_turns = 20;                        % Number of turns p/ magnet [-]
    par.n_pipes = 4;                 % Number of pipes [-]
       
    Cc = 1e-3;      
    w_cb = 10e-3;                   % 
    
    par.NFour = 100;
    tau_c = 4*par.tau_m/3;   
    par.tau = 3*tau_c/2;
    Tau_a = 0.5*tau_c;
    Tau_b = 1.5*tau_c;
    Tau_c = 2.5*tau_c;

    %----------- Fixed relations for Coil offsets
    par.x1a = -tau_c/2 + Cc/2 + Tau_a;
    par.x2a = -tau_c/2 + Cc/2 + w_cb + Tau_a;
    par.x3a = tau_c/2 - Cc/2 - w_cb + Tau_a;
    par.x4a = tau_c/2 - Cc/2 + Tau_a;

    par.x1b = -tau_c/2 + Cc/2 + Tau_b;
    par.x2b = -tau_c/2 + Cc/2 + w_cb + Tau_b;
    par.x3b = tau_c/2 - Cc/2 - w_cb + Tau_b;
    par.x4b = tau_c/2 - Cc/2 + Tau_b;

    par.x1c = -tau_c/2 + Cc/2 + Tau_c;
    par.x2c = -tau_c/2 + Cc/2 + w_cb + Tau_c;
    par.x3c = tau_c/2 - Cc/2 - w_cb + Tau_c;
    par.x4c = tau_c/2 - Cc/2 + Tau_c;
end