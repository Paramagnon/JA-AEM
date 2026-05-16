function Geom = returnMagnetTest(options)
    %--- 4 magnets repeating test
    x0 = 0;
    y0 = 0;
    tau_m = 30e-3;
    w_m = 26e-3;
    h_bi = 5e-3;    
    h_m = 5e-3;
    h_air = 25e-3;
    
    %--- Surface charge
    mu0 = 4*pi*1e-7;
    mur = 1;
    rho = 0e6;             % A/m
    J = 100e6;                % A/m^2
    A = 0;                  % Top boundary condition for testing
    % a = mur*mu0;           % Use for magnets -> solves magnetic potential
    a = 1;                  % Use for coils -> solves magnetic vector potential
    %------------ Domain properties ------------
    Geom(1).a = a; Geom(1).source = -mur*mu0*J;   % Domain property
    Geom(2).a = a; Geom(2).source = mur*mu0*J;   % Domain property
    Geom(3).a = a; Geom(3).source = -mur*mu0*J;   % Domain property
    Geom(4).a = a; Geom(4).source = mur*mu0*J;    % Domain property
    Geom(5).a = a; Geom(5).source = 0;      % Domain property
    
    %------------ Boundaries ------------

    %--- Magnet 1 ---
    p.x0 = x0;
    p.x1 = x0+w_m;
    p.y0 = h_bi;
    p.y1 = h_bi+h_m;
    bndr.idx = [1,3,3,2]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,rho,0];
    p.domain_idx = [0,5,5,5]; %[0,2,3,0];
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,2,1,16];
    [Geom(1).boundary(1)] = createSquareBoundary(p,bndr,options,0);
   
    %--- Magnet 2 ---
    p.x0 = x0+tau_m;
    p.x1 = x0+tau_m+w_m;
    p.y0 = h_bi;
    p.y1 = h_bi+h_m;
    bndr.idx = [1,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr    
    bndr.val = [0,0,-rho,0];
    p.domain_idx = [0,5,5,5]; %[0,2,3,0];
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,6,5,4];
    [Geom(2).boundary(1)] = createSquareBoundary(p,bndr,options,0);

    %--- Magnet 3 ---
    p.x0 = x0+2*tau_m;
    p.x1 = x0+2*tau_m+w_m;
    p.y0 = h_bi;
    p.y1 = h_bi+h_m;
    bndr.idx = [1,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr    
    bndr.val = [0,0,rho,0];
    p.domain_idx = [0,5,5,5]; %[0,2,3,0];
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,10,9,8];
    [Geom(3).boundary(1)] = createSquareBoundary(p,bndr,options,0);

    %--- Magnet 4 ---
    p.x0 = x0+3*tau_m;
    p.x1 = x0+3*tau_m+w_m;
    p.y0 = h_bi;
    p.y1 = h_bi+h_m; 
    bndr.idx = [1,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr    
    bndr.val = [0,0,-rho,0];
    p.domain_idx = [0,5,5,5];
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,14,13,12];    
    [Geom(4).boundary(1)] = createSquareBoundary(p,bndr,options,0);

    %--- Air region ---
    p.x0 = x0;
    p.boundary(1).y0 = y0;
    p.xgeom = [x0,x0+w_m,x0+w_m,x0+tau_m,x0+tau_m,x0+tau_m+w_m,...
        x0+tau_m+w_m,x0+2*tau_m,x0+2*tau_m,x0+2*tau_m+w_m,x0+2*tau_m+w_m,x0+3*tau_m,...
        x0+3*tau_m,x0+3*tau_m+w_m,x0+3*tau_m+w_m,x0+4*tau_m,x0+4*tau_m,x0+4*tau_m,x0,x0];        
    p.ygeom = [y0+h_bi+h_m,y0+h_bi+h_m,y0+h_bi,y0+h_bi,y0+h_bi+h_m,...
        y0+h_bi+h_m,y0+h_bi,y0+h_bi,y0+h_bi+h_m,y0+h_bi+h_m,y0+h_bi,y0+h_bi,y0+h_bi+h_m,...
        y0+h_bi+h_m,y0+h_bi,y0+h_bi,y0+h_bi+h_m,y0+h_bi+h_m+h_air,y0+h_bi+h_m+h_air,y0+h_bi+h_m];

    bndr.idx = [3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,2,2,1,2]; % idx = 1 is u bndr, idx = 0 is un bndr   
    bndr.val = [rho,0,0,0,-rho,0,0,0,rho,0,0,0,-rho,0,0,0,0,A,0];
    p.domain_idx = [1,1,0,2,2,2,0,3,3,3,0,4,4,4,0,1,5,0,5];
    p.boundary_idx = [1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,0,1];
    p.element_idx = [3,2,0,4,3,2,0,4,3,2,0,4,3,2,0,4,19,0,17];   
    [Geom(5).boundary(1)] = createRandomBoundary(p,bndr,options);
end