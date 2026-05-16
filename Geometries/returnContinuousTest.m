function [Geom] = returnContinuousTest(options)
    %--- Continuous boundary test
    x0 = 0;
    y0 = 0;
    w_b = 30e-3;   
    h_b = 30e-3;
  
    %------------ Domain properties ------------
    Geom(1).kappa_f = @kappa_copper; Geom(1).kappa_typ = 50;
    Geom(2).kappa_f = @kappa_steel; Geom(2).kappa_typ = 0; 
    
    Geom(1).sigma_f = @sigma_copper;
    Geom(2).sigma_f = @sigma_steel;

    Geom(1).Cp_f = @Cp_copper;
    Geom(2).Cp_f = @Cp_steel;

    Geom(1).mu_f = @mu_copper;
    Geom(2).mu_f = @mu_steel;
    Geom(2).SCpol = [-1, 1]; 
    %------------ Boundaries ------------

    %--- Square 1 ---
    p.x0 = x0;
    p.x1 = x0+w_b;
    p.y0 = y0;
    p.y1 = y0+h_b;
    bndr.idx = [1,2,1,2]; % > 2 for shared, 4 for SC
    bndr.val = [0,0,0,0];
    bndr.C1 = [1,1,1,-1];
    bndr.C2 = [0,1,0,1]/Geom(1).mu_f(270);
    bndr.C3 = [0,0,0,0];

    bndr.C1T = [1,1,1,-1];
    bndr.C2T = [0,1,0,1]/Geom(1).mu_f(270);
    bndr.C3T = [0,0,0,0];
    p.domain_idx = [0,1,0,1]; 
    p.boundary_idx = [0,1,0,1];
    p.element_idx = [0,4,0,2];
    [Geom(1).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);  
    %--- Inner boundary ---
    p.x0 = x0+0.4*w_b;
    p.x1 = x0+0.6*w_b;
    p.y0 = y0+0.4*h_b;
    p.y1 = y0+0.6*h_b;
    bndr.idx = [3,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,0,0];
    bndr.C1 = [1,1,1,1];
    bndr.C2 = [1,1,1,1]/Geom(1).mu_f(270);
    bndr.C3 = [0,0,0,0];

    bndr.C1T = [1,1,1,1];
    bndr.C2T = [1,1,1,1]/Geom(1).mu_f(270);
    bndr.C3T = [0,0,0,0];
    p.domain_idx = [2,2,2,2]; 
    p.boundary_idx = [1,1,1,1];
    p.element_idx = [4,3,2,1];
    [Geom(1).boundary(2)] = createSquareBoundaryAEM(p,bndr,options,1);

    %--- Square 2 ---
    p.x0 = x0+0.4*w_b;
    p.x1 = x0+0.6*w_b;
    p.y0 = y0+0.4*h_b;
    p.y1 = y0+0.6*h_b;
    bndr.idx = [3,4,3,4]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,0,0];
    bndr.C1 = [-1,-1,-1,-1];
    bndr.C2 = [1,1,1,1]/Geom(2).mu_f(270);
    bndr.C3 = [0,0,0,0];

    bndr.C1T = [-1,-1,-1,-1];
    bndr.C2T = [1,1,1,1]/Geom(2).mu_f(270);
    bndr.C3T = [0,0,0,0];
    p.domain_idx = [1,1,1,1];
    p.boundary_idx = [2,2,2,2];
    p.element_idx = [4,3,2,1];
    [Geom(2).boundary(1)] = createSquareBoundaryAEM(p,bndr,options,0);

    [Geom(1).domain, Geom(1).Ldmn] = insertDomainPoints(Geom(1), options);
    options.dL = 4*options.dL;
    options.rmin = 0.0005;
    [Geom(2).domain, Geom(2).Ldmn] = insertDomainPoints(Geom(2), options);
    
end