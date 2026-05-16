function Geom = returnCoilSolver(options,par)
    %------------ Parameters ------------
    x0 = par.x0;
    y0 = par.y0;
        
    h_c = par.h_c/2;                    % Coil height/2 (symmetry in z-axis) [m]
    h_g = par.h_g;                      % Air gap height [m]
    h_cw = par.h_cw;                    % Cryo wall height [m]
    h_turn = par.h_turn;                % Height per turn [m]
    h_p = par.h_p;                      % Cold plate height [m]

    tau_m = par.tau_m;                  % Magnet period [m]     
    w_turn = par.w_turn;                % Width per turn [m]
    r_pipes = par.r_pipes;              % Radius of pipes [m]

    n_m = par.n_m;                      % Number of magnets [-]
    n_turns = par.n_turns;              % Number of turns p/ magnet [-]
    n_pipes = par.n_pipes;              % Number of pipes [-]

    L = n_m*tau_m;
    dL = L/(n_pipes+1);
   
    dLnodes = options.dL;
    dLnodesx = options.dLx;
    dLnodesy = options.dLy;
    dNnodesy = options.Ny;
    %------------ Boundaries ------------
    %------------ Domain properties ------------
    % Geom(1) - Coil region
    % Geom(2) - Air gap
    % Geom(3) - Cryostat wall
    % Geom(4:4+(Ntapes*2+1)*Nm) - Copper plus air core
    % Geom(4+(Ntapes*2+1)*Nm+1) - Cooling plate 
    % Geom(4+(Ntapes*2+1)*Nm+2:4+(Ntapes*2+1)*Nm+2+Ntubes) - Cooling tubes

    Geom(1).kappa_f = @kappa_copper; Geom(1).kappa_typ = 50;
    Geom(2).kappa_f = @kappa_air; Geom(2).kappa_typ = 0; 
    Geom(3).kappa_f = @kappa_steel; Geom(3).kappa_typ = 0; 

    Geom(1).sigma_f = @sigma_air;
    Geom(2).sigma_f = @sigma_air;
    Geom(3).sigma_f = @sigma_copper;

    Geom(1).Cp_f = @Cp_copper;
    Geom(2).Cp_f = @Cp_air;
    Geom(3).Cp_f = @Cp_steel;

    Geom(1).mu_f = @mu_copper;
    Geom(2).mu_f = @mu_air;
    Geom(3).mu_f = @mu_steel;

    idx = 0;
    for ll = 1:n_m
        % Left HTS tapes
        for nn = 1:n_turns
            Geom(3+idx+nn).kappa_f = @kappa_copper;
            Geom(3+idx+nn).kappa_typ = 50;
            Geom(3+idx+nn).Cp_f = @Cp_copper;
            Geom(3+idx+nn).mu_f = @mu_copper;
            Geom(3+idx+nn).sigma_f = @sigma_copper;
            Geom(3+idx+nn).SCpol = (-1)^mod(ll,2);
        end
        idx = idx + n_turns;
        % Core between tapes
        Geom(3+idx+1).kappa_f = @kappa_air; % --> Replace with epoxy
        Geom(3+idx+1).kappa_typ = 0;
        Geom(3+idx+1).Cp_f = @Cp_air; % --> Replace with epoxy
        Geom(3+idx+1).mu_f = @mu_air; % --> Replace with epoxy
        Geom(3+idx+1).sigma_f = @sigma_air; % --> Replace with epoxy

        idx = idx + 1;
        % Right HTS tapes
        for nn = 1:n_turns
            Geom(3+idx+nn).kappa_f = @kappa_copper;
            Geom(3+idx+nn).kappa_typ = 50;
            Geom(3+idx+nn).Cp_f = @Cp_copper;
            Geom(3+idx+nn).mu_f = @mu_copper;
            Geom(3+idx+nn).sigma_f = @sigma_copper;
            Geom(3+idx+nn).SCpol = (-1)^(mod(ll,2)+1);
        end
        idx = idx + n_turns;
    end

    Geom(3+idx+1).kappa_f = @kappa_copper;
    Geom(3+idx+1).kappa_typ = 50;
    Geom(3+idx+1).Cp_f = @Cp_copper;
    Geom(3+idx+1).mu_f = @mu_copper;
    Geom(3+idx+1).sigma_f = @sigma_copper;

    for nn = 1:n_pipes
        Geom(3+idx+1+nn).kappa_f = @kappa_air;
        Geom(3+idx+1+nn).kappa_typ = 0;
        Geom(3+idx+1+nn).Cp_f = @Cp_air;
        Geom(3+idx+1+nn).mu_f = @mu_air;
        Geom(3+idx+1+nn).sigma_f = @sigma_air;
    end 
    
    %------------------------- Coil region -------------------------
    p.x0 = x0;
    p.x1 = x0+L;
    p.y0 = y0;
    p.y1 = y0+h_c;
    %--- idx 2 for internal continuous bndr, 3 for other geom
    bndr.idx = [0,2,3,2]; 
    bndr.C1 = [0,1,1,-1];
    bndr.C2 = [1,1,1,1]/Geom(1).mu_f(par.T0);
    bndr.C3 = [0,0,0,0];
    % Thermal model BC
    bndr.C1T = [0,1,1,-1];
    bndr.C2T = [1,1,1,1]*Geom(1).kappa_f(par.T0,Geom(1).kappa_typ);
    bndr.C3T = [0,0,0,0]; % ---> check thermal BC

    p.domain_idx = [0,1,2,1];
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,4,1,2];
    [Geom(1).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);
    options.dLx = 0;
    options.dLy = 0;
    options.rmin = max(4/(7*options.Nx),4/(7*options.Ny));
    [Geom(1).domain, Geom(1).Ldmn] = insertDomainPoints(Geom(1), options);
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    %------------------------- Air gap -------------------------
    p.x0 = x0;
    p.x1 = x0+L;
    p.y0 = y0+h_c;
    p.y1 = y0+h_c+h_g;
    %--- idx 2 for internal continuous bndr, 3 for other geom
    bndr.idx = [3,2,3,2]; 
    bndr.C1 = [-1,1,1,-1];
    bndr.C2 = [1,1,1,1]/Geom(2).mu_f(par.T0);
    bndr.C3 = [0,0,0,0];
    % Thermal model BC
    bndr.C1T = [-1,1,1,-1];
    bndr.C2T = [1,1,1,1]*Geom(2).kappa_f(par.T0,Geom(2).kappa_typ);
    bndr.C3T = [0,0,0,0];

    p.domain_idx = [1,2,3,2];
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [3,4,1,2];
    [Geom(2).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);
    options.dLx = 0;
    options.dLy = 0;    
    options.rmin = max(4/(7*options.Nx),4/(7*options.Ny));
    [Geom(2).domain, Geom(2).Ldmn] = insertDomainPoints(Geom(2), options);    
    options.dLx = dLnodesx;
    options.dLy = dLnodesy;
    %----------------------- HTS Tape magnets -----------------------    
    x_hts = [];
    n_dmn = (2*n_turns+1)*n_m;
    options.Ny = 2.5*dNnodesy;
    options.rmin = max(4/(7*options.Nx),4/(7*options.Ny));
    idx0 = 3;
    idx = idx0;
    for ll = 1:n_m
        for nn = 1:n_turns
            idx = idx + 1;

            p.x0 = x0 + (ll-1)*tau_m + (nn-1)*w_turn;
            x_hts(end+1) = p.x0;
            p.x1 = x0 + (ll-1)*tau_m + nn*w_turn;
            p.y0 = y0 + h_c + h_g + h_cw;
            p.y1 = y0 + h_c + h_g + h_cw + h_turn;
            %--- idx 2 for internal continuous bndr, 3 for other geom           
            if (ll == 1 && nn == 1)
                bndr.C1 = [-1,1,1,1];
                bndr.C2 = [1 ,1,1,1]/Geom(idx).mu_f(par.T0);
                bndr.C3 = [0 ,0,0,0];
                % Thermal model BC
                bndr.C1T = [-1,1,1,1];
                bndr.C2T = [1 ,1,1,1]*Geom(idx).kappa_f(par.T0,Geom(idx).kappa_typ);
                bndr.C3T = [0 ,0,0,0];

                p.domain_idx = [idx0,idx+1,idx0+n_dmn+1,idx0+n_dmn];                
                p.element_idx = [n_dmn+2,4,1,2];
            else            
                bndr.C1 = [-1,1,1,-1];
                bndr.C2 = [ 1,1,1, 1]/Geom(idx).mu_f(par.T0);
                bndr.C3 = [ 0,0,0, 0];
                % Thermal model BC
                bndr.C1T = [-1,1,1,-1];
                bndr.C2T = [ 1,1,1, 1]*Geom(idx).kappa_f(par.T0,Geom(idx).kappa_typ);
                bndr.C3T = [ 0,0,0, 0];

                p.domain_idx = [idx0,idx+1,idx0+n_dmn+1,idx-1];
                p.element_idx = [n_dmn+2-(idx-idx0-1),4,idx-idx0,2];
            end
            if mod(nn+1,2)
                bndr.idx = [3,3,3,3+options.SC];
            else
                bndr.idx = [3,3,3,3];
            end
       
            p.boundary_idx = [1,1,1,1];            
            [Geom(idx).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);           
            [Geom(idx).domain, Geom(idx).Ldmn] = insertDomainPoints(Geom(idx), options);
        end
        idx = idx + 1;                

        p.x0 = x0 + (ll-1)*tau_m + n_turns*w_turn;
        x_hts(end+1) = p.x0;
        p.x1 = x0 + ll*tau_m - n_turns*w_turn;
        p.y0 = y0 + h_c + h_g + h_cw;
        p.y1 = y0 + h_c + h_g + h_cw + h_turn;
        %--- idx 2 for internal continuous bndr, 3 for other geom
        bndr.idx = [3,3,3,3];
        bndr.C1 = [-1,1,1,-1];
        bndr.C2 = [ 1,1,1, 1]/Geom(idx).mu_f(par.T0);
        bndr.C3 = [ 0,0,0, 0];
        % Thermal model BC
        bndr.C1T = [-1,1,1,-1];
        bndr.C2T = [ 1,1,1, 1]*Geom(idx).kappa_f(par.T0,Geom(idx).kappa_typ);
        bndr.C3T = [ 0,0,0, 0];

        p.domain_idx = [idx0,idx+1,idx0+n_dmn+1,idx-1];
        p.boundary_idx = [1,1,1,1];
        p.element_idx = [n_dmn+2-(idx-idx0-1),4,idx-idx0,2];
        [Geom(idx).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);
        options.dLx = 0;
        options.dLy = 0;
        [Geom(idx).domain, Geom(idx).Ldmn] = insertDomainPoints(Geom(idx), options);
        options.dLx = dLnodesx;
        options.dLy = dLnodesy;

        for nn = 1:n_turns
            idx = idx + 1;            

            p.x0 = x0 + ll*tau_m - n_turns*w_turn + (nn-1)*w_turn;
            x_hts(end+1) = p.x0;
            p.x1 = x0 + ll*tau_m - n_turns*w_turn + (nn)*w_turn;
            p.y0 = y0 + h_c + h_g + h_cw;
            p.y1 = y0 + h_c + h_g + h_cw + h_turn;

            if (ll == n_m && nn == n_turns)
                bndr.C1 = [-1,-1,1,-1];
                bndr.C2 = [1 , 1,1, 1]/Geom(idx).mu_f(par.T0);
                bndr.C3 = [0 , 0,0, 0];
                % Thermal model BC
                bndr.C1T = [-1,-1,1,-1];
                bndr.C2T = [1 , 1,1, 1]*Geom(idx).kappa_f(par.T0,Geom(idx).kappa_typ);
                bndr.C3T = [0 , 0,0, 0];

                p.domain_idx = [idx0,idx0+1,idx0+n_dmn+1,idx-1];
                p.element_idx = [3,4,n_dmn,2];
            else            
                bndr.C1 = [-1,1,1,-1];
                bndr.C2 = [ 1,1,1, 1]/Geom(idx).mu_f(par.T0);
                bndr.C3 = [ 0,0,0, 0];
                % Thermal model BC
                bndr.C1T = [-1,1,1,-1];
                bndr.C2T = [ 1,1,1, 1]*Geom(idx).kappa_f(par.T0,Geom(idx).kappa_typ);
                bndr.C3T = [ 0,0,0, 0];

                p.domain_idx = [idx0,idx+1,idx0+n_dmn+1,idx-1];
                p.element_idx = [n_dmn+2-(idx-idx0-1),4,idx-idx0,2];
            end
            if mod(nn+1,2)
                bndr.idx = [3,3,3,3+options.SC];
            else
                bndr.idx = [3,3,3,3];
            end
            p.boundary_idx = [1,1,1,1];

            [Geom(idx).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);
            [Geom(idx).domain, Geom(idx).Ldmn] = insertDomainPoints(Geom(idx), options);
        end
    end
    x_hts(end+1) = p.x1;
    idx = idx+1;

    options.Ny = dNnodesy;
    L_hts = length(x_hts);
    %------------------------- Cryo wall -------------------------
    p.x0 = x0;
    p.xgeom = [x_hts(1), x_hts(end), x_hts(end:-1:1), x_hts(1)];
    p.y0 = y0+h_c+h_g;
    p.ygeom = [p.y0,p.y0,ones(1,L_hts)*(p.y0+h_cw),p.y0];
    %--- idx 2 for internal continuous bndr, 3 for other geom
    bndr.idx = [3,2,3*ones(1,L_hts-1),2]; 
    bndr.C1 = [-1,1,ones(1,L_hts-1),-1];
    bndr.C2 = [1,1,ones(1,L_hts-1),1]/Geom(3).mu_f(par.T0);
    bndr.C3 = [0,0,zeros(1,L_hts-1),0];
    % Thermal model BC
    bndr.C1T = [-1,1,ones(1,L_hts-1),-1];
    bndr.C2T = [1,1,ones(1,L_hts-1),1]*Geom(3).kappa_f(par.T0,Geom(3).kappa_typ);
    bndr.C3T = [0,0,zeros(1,L_hts-1),0];

    p.domain_idx = [2,3,L_hts+2:-1:4,3];
    p.boundary_idx = [1,1,ones(1,L_hts-1),1];
    p.element_idx = [3,L_hts+2,ones(1,L_hts-1),2];
    [Geom(3).boundary(1)] = createRandomBoundaryAEM(p,bndr,options);
    options.rmin = max(4/(7*options.Nx),4/(7*options.Ny));
    [Geom(3).domain, Geom(3).Ldmn] = insertDomainPoints(Geom(3), options);

    %----------------------- Cooling pipes -----------------------   
    p.x0 = x0;
    p.xgeom = [x_hts, x_hts(end), x_hts(1), x_hts(1)];
    p.y0 = y0 + h_c + h_g + h_cw + h_turn;
    p.ygeom = [ones(1,L_hts)*p.y0,(p.y0+h_p),(p.y0+h_p),p.y0];   
    %--- idx 2 for internal continuous bndr, 3 for other geom
    bndr.idx = [3*ones(1,L_hts-1),2,0,2];
    bndr.C1 = [-ones(1,L_hts-1),1,1,-1];
    bndr.C2 = [ones(1,L_hts-1),1,0,1]/Geom(idx).mu_f(par.T0);
    bndr.C3 = [zeros(1,L_hts-1),0,0,0];
    % Thermal model BC
    bndr.C1T = [-ones(1,L_hts-1),1,1,-1];
    bndr.C2T = [ones(1,L_hts-1),1,0,1]*Geom(idx).kappa_f(par.T0,Geom(idx).kappa_typ);
    bndr.C3T = [zeros(1,L_hts-1),0,0,0];

    p.domain_idx = [4:n_dmn+3,n_dmn+4,0,n_dmn+4];
    p.boundary_idx = [ones(1,n_dmn),1,0,1];
    p.element_idx = [3*ones(1,n_dmn),n_dmn+3,0,n_dmn+1];
    [Geom(idx).boundary(1)] = createRandomBoundaryAEM(p,bndr,options);

    %--- Cooling pipes ---
    for nn = 1:n_pipes
        p.x0 = x0+nn*dL;
        p.y0 = y0 + h_c + h_g + h_cw + h_turn + h_p/2;
        p.r = r_pipes;
        %--- idx 2 for internal continuous bndr, 3 for other geom
        bndr.idx = 3;
        bndr.C1 = 1;
        bndr.C2 = 1/Geom(idx).mu_f(par.T0);
        bndr.C3 = 0;
        % Thermal model BC
        bndr.C1T = 1;
        bndr.C2T = 1*Geom(idx).kappa_f(par.T0,Geom(idx).kappa_typ);
        bndr.C3T = 0;

        p.domain_idx = idx+nn;
        p.boundary_idx = 1;
        p.element_idx = 1;
        [Geom(idx).boundary(1+nn)] = createCircleBoundaryAEM(p,bndr,options,1);       
    end
    options.rmin = max([4/(7*options.Nx),4/(7*options.Ny),4/(7*options.N)]);
    [Geom(idx).domain, Geom(idx).Ldmn] = insertDomainPoints(Geom(idx), options);

    %--- Add domains within cooling pipes since magnetics continue
    % options.dL = options.dL*3;
    options.dLx = 0;
    options.dLy = 0;
    options.dL = 0;
    for nn = 1:n_pipes
        p.x0 = x0+nn*dL;
        p.y0 = y0 + h_c + h_g + h_cw + h_turn + h_p/2;
        p.r = r_pipes;
        %--- idx 2 for internal continuous bndr, 3 for other geom
        bndr.idx = 3;
        bndr.C1 = -1;
        bndr.C2 = 1/Geom(idx+nn).mu_f(par.T0);
        bndr.C3 = 0;
        % Thermal model BC
        bndr.C1T = -1;
        bndr.C2T = 1*Geom(idx+nn).kappa_f(par.T0,Geom(idx+nn).kappa_typ);
        bndr.C3T = 0;

        p.domain_idx = idx;
        p.boundary_idx = 1+nn;
        p.element_idx = 1;
        [Geom(idx+nn).boundary(1)] = createCircleBoundaryAEM(p,bndr,options,0);       
        options.rmin = 4/(7*options.N);
        [Geom(idx+nn).domain, Geom(idx+nn).Ldmn] = insertDomainPoints(Geom(idx+nn), options);
    end
end