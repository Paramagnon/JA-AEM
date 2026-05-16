function par = loadParCoilCOMSOLYoke()
    par.x0 = 0;
    par.y0 = 0;
        
    % wc	alphac*tau_c	0.012 m	Coilside width
    % id	tau_c-dc-2*wc	0.011 m	Inner diameter of the coil
    % dc	5 [mm]	0.005 m	Distance between consecutive coilsides
    % hc	10 [mm]	0.01 m	Coil height
    % hce2	2 [mm]	0.002 m	Plate thickness active shield
    % hce	4 [mm]	0.004 m	Conductive layer height (Cryostat wall)
    % D	140 [mm]	0.14 m	Depth of the domain (z-direction)
    % gl_m	4 [mm]	0.004 m	Mechanical airgap
    % gl_c	1.5 [mm] - (hce2 -2 [mm])/2	0.0015 m	Airgap either side active shield
    % tau_c	4*tau_m/3	0.04 m	Coil Period
    % h_inf	3 [mm]	0.003 m	Height of infinite regions
    % h_cp	4 [mm]	0.004 m	Height cold plate
    % g	2*gl_c+gl_m+hce+hce2	0.013 m
    % alphac	0.3	0.3
    % fill	1	1	Fill factor of coil
    % sp_turn	0.5 [mm]	5E-4 m	Separation Between Tapes
    % dsc_coil	nr_turns*sp_turn/fill	0.0025 m	
    % dsc	0.002 [mm]	2E-6 m	Tape thickness 0.002 [mm]
    % hsc	8 [mm]	0.008 m	Tape Width (Height = y)
    % id_sc	20 [mm]	0.02 m	SC Coil Inner Diameter
    % xip_sc	(tau_m-id_sc-2*dsc_coil)/2	0.0025 m	Initial Position of the SC Coil (dc/3 )
    % xin_sc	xip_sc+nr_turns*sp_turn/fill+id_sc	0.025 m	
    % tau_m	0.03 [m]	0.03 m	
    % nr_turns	5	5	Turn Number of the SC Coil

    par.h_c = 10e-3;                    % Coil height/2 (symmetry in z-axis) [m]
    par.h_g = 1e-3;                      % Air gap height [m]
    par.h_cw = 5e-3;                    % Cryo wall height [m]
    par.h_turn = 4e-3;                      % Height per turn [m]
    par.h_p = 5e-3;                
    
    par.tau_m =  40e-3;                  % Magnet period [m] 
    par.w_turn = 114e-6;                    % Width per turn [m]
    
    par.r_pipes = 1e-3;                 % Radius of pipes [m]

    par.n_m = 4;                            % Number of magnets [-]
    par.n_turns = 2;                        % Number of turns p/ magnet [-]
    par.n_pipes = 4;                 % Number of pipes [-]
       
    Cc = 5e-3;      
    w_cb = 10e-3;                   % 
    i_d = 20e-3;

    par.NFour = 150;
    tau_c = Cc + 2*w_cb + i_d;
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