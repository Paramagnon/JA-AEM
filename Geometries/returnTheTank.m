function Geom = returnTheTank(options,par)
    %--- Cold plate parameters
    x0 = par.x0;
    y0 = par.y0;
    L = par.L;             % Length in [m]
    H = par.H;              % Height in [m]    
    Npipe = par.Npipe;              % Number of pipes

    T1 = par.T1;
    T2 = par.T2;
    r = par.r;
    %--- Properties
    a = 1;
    %------------ Domain properties ------------
    Geom(1).a = a; Geom(1).source = 0;   % Domain property
    %------------ Boundaries ------------

    %--- Outer Boundary ---
    p.x0 = x0;
    p.x1 = x0+L;
    p.y0 = y0;
    p.y1 = y0+H;
    bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [T1,T2,T2,T2];
    p.domain_idx = [0,0,0,0]; 
    p.boundary_idx = [0,0,0,0];
    p.element_idx = [0,0,0,0];
    [Geom(1).boundary(1)] = createSquareBoundary(p,bndr,options,0);

    %--- Cooling pipes --- 
    dL = L/(Npipe+1);
    for nn = 1:Npipe
        p.x0 = x0+nn*dL;    
        p.y0 = y0+H/2; 
        p.r = r;
        bndr.idx = [1]; % idx = 1 is u bndr, idx = 0 is un bndr    
        bndr.val = [100];
        p.domain_idx = [0]; %[0,2,3,0];
        p.boundary_idx = [0];
        p.element_idx = [0];
        [Geom(1).boundary(1+nn)] = createCircleBoundary(p,bndr,options,1);  

    end



    % Nm = par.Nm;
    % Nturns = par.Nturns;
    % Taum = par.Taum;
    % wm = par.wm;

    Nm = 4;
    Nturns = 5;
    Taum = L/Nm;
    wturn = 140e-6; 
    hturn = 4e-3;
    options.N = options.N*10;
    idx = 1;
    for ll = 1:Nm
        for nn = 1:Nturns
            idx = idx + 1;

            p.x0 = x0 + (ll-1)*Taum + (nn-1)*wturn;
            p.x1 = x0 + (ll-1)*Taum + nn*wturn;
            p.y0 = y0 + H;
            p.y1 = y0 + H + hturn;

            bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
            bndr.val = [T1,T2,T2,T2];
            p.domain_idx = [0,0,0,0]; 
            p.boundary_idx = [0,0,0,0];
            p.element_idx = [0,0,0,0];
            [Geom(idx).boundary(1)] = createSquareBoundary(p,bndr,options,0);
        end

        idx = idx + 1;

        p.x0 = x0 + (ll-1)*Taum + Nturns*wturn;
        p.x1 = x0 + ll*Taum - Nturns*wturn;
        p.y0 = y0 + H;
        p.y1 = y0 + H + hturn;

        bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
        bndr.val = [T1,T2,T2,T2];
        p.domain_idx = [0,0,0,0];
        p.boundary_idx = [0,0,0,0];
        p.element_idx = [0,0,0,0];
        [Geom(idx).boundary(1)] = createSquareBoundary(p,bndr,options,0);

        for nn = 1:Nturns
            idx = idx + 1;

            p.x0 = x0 + ll*Taum - Nturns*wturn + (nn-1)*wturn;
            p.x1 = x0 + ll*Taum - Nturns*wturn + (nn)*wturn;
            p.y0 = y0 + H;
            p.y1 = y0 + H + hturn;

            bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
            bndr.val = [T1,T2,T2,T2];
            p.domain_idx = [0,0,0,0]; 
            p.boundary_idx = [0,0,0,0];
            p.element_idx = [0,0,0,0];
            [Geom(idx).boundary(1)] = createSquareBoundary(p,bndr,options,0);
        end

    end
















end