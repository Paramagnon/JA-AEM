function Geom = returnCoilCOMSOLYoke(options,par)
    %------------ Parameters ------------
    x0 = 0;
    y0 = 0;
        
    h_inf = 20*3e-3;
    h_air_top = 4e-3;
    h_ab = 3.5e-3;
    h_c = 10e-3;                        % Coil height/2 (symmetry in z-axis) [m]
    h_g = 2.5e-3+2.5e-3;                % Air gap height [m]
    h_cw = 3.5e-3;                      % Cryo wall height [m]
    h_turn = 4e-3;                     % Height per turn [m]   

    w_turn = 0.6e-3;                % Width per turn [m]    
    w_left = 0.05;
    w_mid = 2.5e-3; 

    tau = 67.5e-3;
    L = 2*tau; 
   
    fracrmin = 1; % 4/7 is original

    dLnodesx = options.dLx;
    dLnodesy = options.dLy;
    dNnodesy = options.Ny;
    %------------ Domain properties ------------
    % Geom(1) - Air region h_inf + h_ab
    % Geom(2) - Coil region
    % Geom(3) - Air gap
    % Geom(4) - Conductive wall
    % Geom(5) - Air gap
    % Geom(4:4+(Ntapes*2+1)*Nm) - Copper plus air core
    % Geom(4:4+(Ntapes*2+1)*Nm+1) - Air top

    Geom(1).kappa_f = @kappa_air; Geom(1).kappa_typ = 0;
    Geom(2).kappa_f = @kappa_air; Geom(2).kappa_typ = 0;
    Geom(3).kappa_f = @kappa_air; Geom(3).kappa_typ = 0; 
    Geom(4).kappa_f = @kappa_copper; Geom(4).kappa_typ = 50; 
    Geom(5).kappa_f = @kappa_copper; Geom(5).kappa_typ = 50; 
    Geom(6).kappa_f = @kappa_copper; Geom(6).kappa_typ = 50; 
    
    Geom(1).sigma_f = @sigma_air;
    Geom(2).sigma_f = @sigma_air;
    Geom(3).sigma_f = @sigma_air;
    Geom(4).sigma_f = @sigma_copper;
    Geom(5).sigma_f = @sigma_air;
    Geom(6).sigma_f = @sigma_air;

    Geom(1).Cp_f = @Cp_air;
    Geom(2).Cp_f = @Cp_air;
    Geom(3).Cp_f = @Cp_air;
    Geom(4).Cp_f = @Cp_copper;
    Geom(5).Cp_f = @Cp_copper;
    Geom(6).Cp_f = @Cp_copper;

    Geom(1).rho_f = @rho_copper;
    Geom(2).rho_f = @rho_copper;
    Geom(3).rho_f = @rho_copper;
    Geom(4).rho_f = @rho_copper;
    Geom(5).rho_f = @rho_copper;
    Geom(6).rho_f = @rho_copper;

    Geom(1).mu_f = @mu_air;
    Geom(2).mu_f = @mu_air;
    Geom(3).mu_f = @mu_air;
    Geom(4).mu_f = @mu_copper;
    Geom(5).mu_f = @mu_air;
    Geom(6).mu_f = @mu_air;
    
    for nn = 1:4
        Geom(6+nn).kappa_f = @kappa_copper;
        Geom(6+nn).kappa_typ = 50;
        Geom(6+nn).Cp_f = @Cp_copper;
        Geom(6+nn).rho_f = @rho_copper;
        Geom(6+nn).mu_f = @mu_air;
        Geom(6+nn).sigma_f = @sigma_air;
        Geom(6+nn).SCpol = 1;
    end
    % Core between tapes
    Geom(11).kappa_f = @kappa_copper; % --> Replace with epoxy
    Geom(11).kappa_typ = 50;
    Geom(11).Cp_f = @Cp_copper; % --> Replace with epoxy
    Geom(11).rho_f = @rho_copper; % --> Replace with epoxy
    Geom(11).mu_f = @mu_air; % --> Replace with epoxy
    Geom(11).sigma_f = @sigma_air; % --> Replace with epoxy
    Geom(11).SCpol = 1;

    % Right HTS tapes
    for nn = 1:4
        Geom(11+nn).kappa_f = @kappa_copper;
        Geom(11+nn).kappa_typ = 50;
        Geom(11+nn).Cp_f = @Cp_copper;
        Geom(11+nn).rho_f = @rho_copper;
        Geom(11+nn).mu_f = @mu_air;
        Geom(11+nn).sigma_f = @sigma_air;
        Geom(11+nn).SCpol = -1;
    end    
        
    % Right block
    Geom(16).kappa_f = @kappa_copper; % --> Replace with epoxy
    Geom(16).kappa_typ = 50;
    Geom(16).Cp_f = @Cp_copper; % --> Replace with epoxy
    Geom(16).rho_f = @rho_copper; % --> Replace with epoxy
    Geom(16).mu_f = @mu_air; % --> Replace with epoxy
    Geom(16).sigma_f = @sigma_air; % --> Replace with epoxy
    Geom(16).SCpol = -1;

    % Top air block
    Geom(17).kappa_f = @kappa_air; % --> Replace with epoxy
    Geom(17).kappa_typ = 0;
    Geom(17).Cp_f = @Cp_air; % --> Replace with epoxy
    Geom(17).rho_f = @rho_copper; % --> Replace with epoxy
    Geom(17).mu_f = @mu_air; % --> Replace with epoxy
    Geom(17).sigma_f = @sigma_air; % --> Replace with epoxy  

    
    %---------------------- Bottom air region ----------------------
    p.x0 = x0;
    p.x1 = x0+L;
    p.y0 = y0-h_ab-h_inf;
    p.y1 = y0;
    % Magnetic model BC
    bndr.idx = [1,2,3,2]; 
    bndr.C1 = [1,1,1,-1];
    bndr.C2 = [0,1,1,1]/Geom(1).mu_f(par.T0);
    bndr.C3 = [0,0,0,0];
    % Thermal model BC
    bndr.idxT = [1,2,3,2]; 
    bndr.C1T = [0,1,1,-1];
    bndr.C2T = [1,1,1,1]*Geom(1).kappa_f(par.T0,Geom(1).kappa_typ);
    bndr.C3T = [0,0,0,0]; % ---> check thermal BC

    p.domain_idx = [0,1,2,1];
    p.domain_idxT = [0,1,2,1];
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,4,1,2];
    options.Ny = dNnodesy/5;
    [Geom(1).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);
    options.Ny = dNnodesy;
    options.dLx = 0;
    options.dLy = 0;
    options.rmin = max(fracrmin/(options.Nx),fracrmin/(options.Ny));
    [Geom(1).domain, Geom(1).Ldmn] = insertDomainPoints(Geom(1), options);
    [Geom(1).domainT, Geom(1).LdmnT] = insertDomainPointsThermal(Geom(1), options, 0);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;

    %------------------------- Coil region -------------------------
    p.x0 = x0;
    p.x1 = x0+L;
    p.y0 = y0;
    p.y1 = y0+h_c;
    % Magnetic model BC
    bndr.idx = [3,2,3,2]; 
    bndr.C1 = [-1,1,1,-1];
    bndr.C2 = [1,1,1,1]/Geom(2).mu_f(par.T0);
    bndr.C3 = [0,0,0,0];
    % Thermal model BC
    bndr.idxT = [3,2,3,2]; 
    bndr.C1T = [0,1,1,-1];
    bndr.C2T = [1,1,1,1]*Geom(2).kappa_f(par.T0,Geom(2).kappa_typ);
    bndr.C3T = [0,0,0,0]; % ---> check thermal BC

    p.domain_idx = [1,2,3,2];
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [3,4,1,2];
    options.dLx = 0;
    options.dLy = 0;
    [Geom(2).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);
    options.rmin = max(fracrmin/(options.Nx),fracrmin/(options.Ny));
    [Geom(2).domain, Geom(2).Ldmn] = insertDomainPoints(Geom(2), options);
    [Geom(2).domainT, Geom(2).LdmnT] = insertDomainPointsThermal(Geom(2), options, 0);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    %------------------------- Air gap -------------------------
    p.x0 = x0;
    p.x1 = x0+L;
    p.y0 = y0+h_c;
    p.y1 = y0+h_c+h_g;
    % Magnetic model BC
    bndr.idx = [3,2,3,2]; 
    bndr.C1 = [-1,1,1,-1];
    bndr.C2 = [1,1,1,1]/Geom(3).mu_f(par.T0);
    bndr.C3 = [0,0,0,0];
    % Thermal model BC
    bndr.idxT = [3,2,3,2]; 
    bndr.C1T = [-1,1,1,-1];
    bndr.C2T = [1,1,1,1]*Geom(3).kappa_f(par.T0,Geom(3).kappa_typ);
    bndr.C3T = [0,0,0,0];

    p.domain_idx = [2,3,4,3];
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [3,4,1,2];
    [Geom(3).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);
    options.dLx = 0;
    options.dLy = 0;    
    options.rmin = max(fracrmin/(options.Nx),fracrmin/(options.Ny));
    [Geom(3).domain, Geom(3).Ldmn] = insertDomainPoints(Geom(3), options);    
    [Geom(3).domainT, Geom(3).LdmnT] = insertDomainPointsThermal(Geom(3), options, 0);    
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;

    %------------------------- Conductive region -------------------------
    p.x0 = x0;
    p.x1 = x0+L;
    p.y0 = y0+h_c+h_g;
    p.y1 = y0+h_c+h_g+h_cw;
    % Magnetic model BC
    bndr.idx = [3,2,3,2]; 
    bndr.C1 = [-1,1,1,-1];
    bndr.C2 = [1,1,1,1]/Geom(4).mu_f(par.T0);
    bndr.C3 = [0,0,0,0];
    % Thermal model BC
    bndr.idxT = [0,2,3,2]; 
    bndr.C1T = [1,1,1,-1];
    bndr.C2T = [0,1,1,1]*Geom(4).kappa_f(par.T0,Geom(4).kappa_typ);
    bndr.C3T = [par.T0,0,0,0]; % ---> check thermal BC
    % bndr.idxT = [1,2,3,2]; 
    % bndr.C1T = [0,1,1,-1];
    % bndr.C2T = [1,1,1,1]*Geom(4).kappa_f(par.T0,Geom(4).kappa_typ);
    % bndr.C3T = [0,0,0,0]; % ---> check thermal BC

    p.domain_idx = [3,4,5,4];
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [3,4,1,2];
    [Geom(4).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);    
    options.rmin = max(fracrmin/(options.Nx),fracrmin/(options.Ny));
    % options.dLx = 0;
    % options.dLy = 0;
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    [Geom(4).domain, Geom(4).Ldmn] = insertDomainPoints(Geom(4), options);        
    [Geom(4).domainT, Geom(4).LdmnT] = insertDomainPointsThermal(Geom(4), options, 0);    

    %----------------------- HTS Tape magnets -----------------------    
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    options.Ny = options.NySC;
    options.rmin = max(fracrmin/(options.Nx),fracrmin/(options.Ny));
    %--- Left air block
    p.x0 = x0;
    p.x1 = x0+w_left;
    p.y0 = y0+h_c+h_g+h_cw+h_ab;
    p.y1 = y0+h_c+h_g+h_cw+h_ab+h_turn;
    % Magnetic model BC
    bndr.idx = [3,3,3,3]; 
    bndr.C1 = [-1,1,1,1];
    bndr.C2 = [1,1,1,1]/Geom(6).mu_f(par.T0);
    bndr.C3 = [0,0,0,0];
    % Thermal model BC
    bndr.idxT = [3,3,0,3]; 
    bndr.C1T = [-1,1,1,1];
    bndr.C2T = [1,1,0,1]*Geom(6).kappa_f(par.T0,Geom(6).kappa_typ);
    bndr.C3T = [0,0,par.T0,0]; % ---> check thermal BC
    % bndr.idxT = [3,3,0,3];
    % bndr.C1T = [-1,1,1,1];
    % bndr.C2T = [1,1,0,1]*Geom(6).kappa_f(par.T0,Geom(6).kappa_typ);
    % bndr.C3T = [0,0,par.T0,0]; % ---> check thermal BC
    options.dLx = 0;
    options.dLy = 0;
    p.domain_idx = [5,7,17,16];
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [13,4,1,2];
    [Geom(6).boundary(1)] = createSquareBoundarySC(p,bndr,options,0);    
    options.rmin = max(fracrmin/(options.Nx),fracrmin/(options.Ny));
    [Geom(6).domain, Geom(6).Ldmn] = insertDomainPoints(Geom(6), options);  
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    [Geom(6).domainT, Geom(6).LdmnT] = insertDomainPointsThermal(Geom(6), options,1);    
    options.dLx = 0;
    options.dLy = 0;

    x_hts = [x0];
    for nn = 1:4        
        p.x0 = x0 + w_left + (nn-1)*w_turn;
        x_hts(end+1) = p.x0;
        p.x1 = x0 + w_left + nn*w_turn;
        p.y0 = y0 + h_c + h_g + h_cw + h_ab;
        p.y1 = y0 + h_c + h_g + h_cw + h_ab + h_turn;
        % Magnetic model BC
        bndr.C1 = [-1,1,1,-1];
        bndr.C2 = [ 1,1,1, 1]/Geom(6+nn).mu_f(par.T0);
        bndr.C3 = [ 0,0,0, 0];
        % Thermal model BC
        bndr.C1T = [-1,1,1,-1];
        bndr.C2T = [ 1,1,0, 1]*Geom(6+nn).kappa_f(par.T0,Geom(6+nn).kappa_typ);
        bndr.C3T = [ 0,0,par.T0, 0];

        p.domain_idx = [5,7+nn,17,5+nn];
        p.element_idx = [13-nn,4,1+nn,2];      
        bndr.idx = [3,3,3,4];   
        bndr.idxT = [3,3,0,3]; 
        p.boundary_idx = [1,1,1,1];
        [Geom(6+nn).boundary(1)] = createSquareBoundarySC(p,bndr,options,0);
        [Geom(6+nn).domain, Geom(6+nn).Ldmn] = insertDomainPoints(Geom(6+nn), options);
        options.dLx = dLnodesx;
        options.dLy = dLnodesy;
        [Geom(6+nn).domainT, Geom(6+nn).LdmnT] = insertDomainPointsThermal(Geom(6+nn), options,1);
        options.dLx = 0;
        options.dLy = 0;
    end    

    p.x0 = x0 + w_left + 4*w_turn;
    x_hts(end+1) = p.x0;
    p.x1 = x0 + w_left + 4*w_turn + w_mid;
    p.y0 = y0 + h_c + h_g + h_cw + h_ab;
    p.y1 = y0 + h_c + h_g + h_cw + h_ab + h_turn;
    % Magnetic model BC
    bndr.idx = [3,3,3,4];
    bndr.C1 = [-1,1,1,-1];
    bndr.C2 = [ 1,1,1, 1]/Geom(11).mu_f(par.T0);
    bndr.C3 = [ 0,0,0, 0];
    % Thermal model BC
    bndr.idxT = [3,3,0,3];
    bndr.C1T = [-1,1,1,-1];
    bndr.C2T = [ 1,1,0, 1]*Geom(11).kappa_f(par.T0,Geom(11).kappa_typ);
    bndr.C3T = [ 0,0,par.T0, 0];

    p.domain_idx = [5,12,17,10];
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [8,4,6,2];
    [Geom(11).boundary(1)] = createSquareBoundarySC(p,bndr,options,0);
    [Geom(11).domain, Geom(11).Ldmn] = insertDomainPoints(Geom(11), options);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    [Geom(11).domainT, Geom(11).LdmnT] = insertDomainPointsThermal(Geom(11), options, 1);
    options.dLx = 0;
    options.dLy = 0;

    for nn = 1:4        
        p.x0 = x0 + w_left + (4+nn-1)*w_turn + w_mid;
        x_hts(end+1) = p.x0;
        p.x1 = x0 + w_left + (4+nn)*w_turn + w_mid;
        p.y0 = y0 + h_c + h_g + h_cw + h_ab;
        p.y1 = y0 + h_c + h_g + h_cw + h_ab + h_turn;
        % Magnetic model BC
        bndr.C1 = [-1,1,1,-1];
        bndr.C2 = [ 1,1,1, 1]/Geom(11+nn).mu_f(par.T0);
        bndr.C3 = [ 0,0,0, 0];
        % Thermal model BC
        bndr.C1T = [-1,1,1,-1];
        bndr.C2T = [ 1,1,0, 1]*Geom(11+nn).kappa_f(par.T0,Geom(6+nn).kappa_typ);
        bndr.C3T = [ 0,0,par.T0, 0];

        p.domain_idx = [5,12+nn,17,10+nn];
        p.element_idx = [8-nn,4,6+nn,2];      
        bndr.idx = [3,3,3,4];    
        bndr.idxT = [3,3,0,3]; 
        p.boundary_idx = [1,1,1,1];
        [Geom(11+nn).boundary(1)] = createSquareBoundarySC(p,bndr,options,0);
        [Geom(11+nn).domain, Geom(11+nn).Ldmn] = insertDomainPoints(Geom(11+nn), options);
        options.dLx = dLnodesx;
        options.dLy = dLnodesy;
        [Geom(11+nn).domainT, Geom(11+nn).LdmnT] = insertDomainPointsThermal(Geom(11+nn), options, 1);
        options.dLx = 0;
        options.dLy = 0;
    end  

    %--- Right air block
    p.x0 = x0 + w_left + 8*w_turn + w_mid;
    x_hts(end+1) = p.x0;
    p.x1 = x0 + L;
    p.y0 = y0 + h_c + h_g + h_cw + h_ab;
    p.y1 = y0 + h_c + h_g + h_cw + h_ab + h_turn;
    % Magnetic model BC
    bndr.idx = [3,3,3,4];
    bndr.C1 = [-1,-1,1,-1];
    bndr.C2 = [ 1,1,1, 1]/Geom(16).mu_f(par.T0);
    bndr.C3 = [ 0,0,0, 0];
    % Thermal model BC
    bndr.idxT = [3,3,0,3];
    bndr.C1T = [-1,-1,1,-1];
    bndr.C2T = [ 1,1,0, 1]*Geom(16).kappa_f(par.T0,Geom(16).kappa_typ);
    bndr.C3T = [ 0,0,par.T0, 0];
    % bndr.idxT = [3,3,0,3];
    % bndr.C1T = [-1,-1,1,-1];
    % bndr.C2T = [ 1,1,0, 1]*Geom(16).kappa_f(par.T0,Geom(16).kappa_typ);
    % bndr.C3T = [ 0,0,par.T0, 0];

    p.domain_idx = [5,6,17,15];
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [3,4,11,2];
    [Geom(16).boundary(1)] = createSquareBoundarySC (p,bndr,options,0);
    [Geom(16).domain, Geom(16).Ldmn] = insertDomainPoints(Geom(16), options);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    [Geom(16).domainT, Geom(16).LdmnT] = insertDomainPointsThermal(Geom(16), options, 1);

    x_hts(end+1) = p.x1;

    options.Ny = dNnodesy;
    L_hts = length(x_hts);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    %----------------------- Air gap Cryo ----------------------
    options.dLx = 0;
    options.dLy = 0;
    p.x0 = x0;
    p.xgeom = [x_hts(1), x_hts(end), x_hts(end:-1:1), x_hts(1)];
    p.y0 = y0+h_c+h_g+h_cw;
    p.ygeom = [p.y0,p.y0,ones(1,L_hts)*(p.y0+h_ab),p.y0];
    % Magnetic model BC
    bndr.idx = [3,2,3*ones(1,L_hts-1),2]; 
    bndr.C1 = [-1,1,ones(1,L_hts-1),-1];
    bndr.C2 = [1,1,ones(1,L_hts-1),1]/Geom(5).mu_f(par.T0);
    bndr.C3 = [0,0,zeros(1,L_hts-1),0];
    % Thermal model BC
    bndr.idxT = [3,2,3*ones(1,L_hts-1),2]; 
    bndr.C1T = [-1,1,ones(1,L_hts-1),-1];
    bndr.C2T = [1,1,ones(1,L_hts-1),1]*Geom(5).kappa_f(par.T0,Geom(5).kappa_typ);
    bndr.C3T = [0,0,zeros(1,L_hts-1),0];
    % bndr.idxT = [3,2,3*ones(1,L_hts-1),2]; 
    % bndr.C1T = [-1,1,ones(1,L_hts-1),-1];
    % bndr.C2T = [1,1,ones(1,L_hts-1),1]*Geom(5).kappa_f(par.T0,Geom(5).kappa_typ);
    % bndr.C3T = [0,0,zeros(1,L_hts-1),0];

    p.domain_idx = [4,5,L_hts+4:-1:6,5];
    p.boundary_idx = [1,1,ones(1,L_hts-1),1];
    p.element_idx = [3,L_hts+2,ones(1,L_hts-1),2];
    [Geom(5).boundary(1)] = createRandomBoundaryAEM(p,bndr,options);
    options.rmin = max(fracrmin/(options.Nx),fracrmin/(options.Ny));
    [Geom(5).domain, Geom(5).Ldmn] = insertDomainPoints(Geom(5), options);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    [Geom(5).domainT, Geom(5).LdmnT] = insertDomainPointsThermal(Geom(5), options, 0);    
    %----------------------- Air Gap top -----------------------   
    options.dLx = 0;
    options.dLy = 0;
    p.x0 = x0;
    p.xgeom = [x_hts, x_hts(end), x_hts(1), x_hts(1)];
    p.y0 = y0 + h_c + h_g + h_cw + h_ab + h_turn;
    p.ygeom = [ones(1,L_hts)*p.y0,(p.y0+h_inf+h_air_top),(p.y0+h_inf+h_air_top),p.y0];   
    % Magnetic model BC
    bndr.idx = [3*ones(1,L_hts-1),2,0,2];
    bndr.C1 = [-ones(1,L_hts-1),1,1,-1];
    bndr.C2 = [ones(1,L_hts-1),1,0,1]/Geom(17).mu_f(par.T0);
    bndr.C3 = [zeros(1,L_hts-1),0,0,0];
    % Thermal model BC
    bndr.idxT = [3*ones(1,L_hts-1),2,0,2];
    bndr.C1T = [-ones(1,L_hts-1),1,1,-1];
    bndr.C2T = [ones(1,L_hts-1),1,0,1]*Geom(17).kappa_f(par.T0,Geom(17).kappa_typ);
    bndr.C3T = [zeros(1,L_hts-1),0,0,0];

    p.domain_idx = [6:16,17,0,17];
    p.boundary_idx = [ones(1,L_hts-1),1,1,1];
    p.element_idx = [3*ones(1,L_hts-1),14,0,12];
    options.Ny = dNnodesy/5;
    [Geom(17).boundary(1)] = createRandomBoundaryAEM(p,bndr,options);
    options.Ny = dNnodesy;
    options.rmin = max([fracrmin/(options.Nx),fracrmin/(options.Ny),fracrmin/(options.N)]);
    [Geom(17).domain, Geom(17).Ldmn] = insertDomainPoints(Geom(17), options);
    [Geom(17).domainT, Geom(17).LdmnT] = insertDomainPointsThermal(Geom(17), options,0);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;


    [Geom] = linkNghbrBndrs(Geom);
    [Geom] = returnNodesAEM(Geom);
    Geom(1).Tlow = 4; % Thermal geom low idx (first thermal)
    Geom(1).Tup = 16; % Thermal geom up idx (last thermal)

end