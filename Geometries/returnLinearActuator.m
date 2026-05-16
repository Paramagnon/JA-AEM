function Geom = returnLinearActuator(options,par)
    %======================  Actuator properties  =========================
    x0 = par.x0;
    y0 = par.y0;
    dx = par.dx; 

    tau_m = par.tau_m; 
    tau_c = par.tau_c;    
    w_m = par.w_m;    
    w_c = par.w_c;
    h_c = par.h_c;
    h_m = par.h_m;   
    h_bi = par.h_bi;
    g = par.g;
    dwm = par.dwm;
    
    if dwm <= 0 
        error("Magnet width can not be equal or larger than magnet period");
    end

    mu0 = 4*pi*1e-7;
    mur = 1.05;    
    %=======================  Domain properties  ==========================
    Geom(1).a = 1.05;   Geom(1).source = -mu0*mur*par.Ju*1e6;   % Domain property
    Geom(2).a = 1.05;   Geom(2).source = -mu0*mur*par.Jv*1e6;   % Domain property
    Geom(3).a = 1.05;   Geom(3).source = -mu0*mur*par.Jw*1e6;   % Domain property
    Geom(4).a = 1;      Geom(4).source = 0;         % Domain property
    Geom(5).a = 1;      Geom(5).source = 0;         % Domain property
    Geom(6).a = 1;      Geom(6).source = 0;         % Domain property
    Geom(7).a = 1;      Geom(7).source = 0;         % Domain property
    Geom(8).a = 1;      Geom(8).source = 0;         % Domain property
    Geom(9).a = 1;      Geom(9).source = 0;         % Domain property
    Geom(10).a = 1;     Geom(10).source = 0;        % Domain property
    Geom(11).a = 1;   Geom(11).source = 0;        % Domain property
    %=============================  Coils  ================================
    %--- Phase A --- 
    p.x0 = x0+tau_m/2+dx;
    p.x1 = x0+tau_m/2+w_c+dx;
    p.y0 = y0;
    p.y1 = y0+h_c/2;
    bndr.idx = [0,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,0,0];
    p.domain_idx = [0,4,4,4]; 
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,4,3,2];
    [Geom(1).boundary(1)] = createSquareBoundary(p,bndr,options,0);          
    
    %--- Phase B --- 
    p.x0 = x0+tau_m/2+tau_c+dx;
    p.x1 = x0+tau_m/2+tau_c+w_c+dx;
    p.y0 = y0;
    p.y1 = y0+h_c/2;
    bndr.idx = [0,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,0,0];
    p.domain_idx = [0,4,4,4]; 
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,8,7,6];
    [Geom(2).boundary(1)] = createSquareBoundary(p,bndr,options,0);          

    %--- Phase C --- 
    p.x0 = x0+tau_m/2+2*tau_c+dx;
    p.x1 = x0+tau_m/2+2*tau_c+w_c+dx;
    p.y0 = y0;
    p.y1 = y0+h_c/2;
    bndr.idx = [0,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,0,0];
    p.domain_idx = [0,4,4,4]; 
    p.boundary_idx = [0,1,1,1];
    p.element_idx = [0,12,11,10];
    [Geom(3).boundary(1)] = createSquareBoundary(p,bndr,options,0);          
    
    %===========================  Air region  =============================
    p.x0 = x0;
    p.y0 = y0;  

    p.xgeom = zeros(1,4+4*6+4*3+1);
    p.ygeom = zeros(1,4+4*6+4*3+1);
    p.domain_idx = zeros(1,4+4*6+4*3);
    p.boundary_idx = zeros(1,4+4*6+4*3);
    p.element_idx = zeros(1,4+4*6+4*3);
    bndr.idx = zeros(1,4+4*6+4*3);
    bndr.val = zeros(1,4+4*6+4*3);       

    p.xgeom(1) = x0;
    p.ygeom(1) = y0;    
    for ii = 1:3        
        p.xgeom(2+(ii-1)*4:2+(ii-1)*4+3) = [x0+tau_m/2+dx+(ii-1)*tau_c,x0+tau_m/2+dx+(ii-1)*tau_c,...
            x0+tau_m/2+dx+(ii-1)*tau_c+w_c,x0+tau_m/2+dx+(ii-1)*tau_c+w_c];
        p.ygeom(2+(ii-1)*4:2+(ii-1)*4+3) = [y0,y0+h_c/2,y0+h_c/2,y0];        
        bndr.idx(2+(ii-1)*4:2+(ii-1)*4+2) = [3,3,3];
        bndr.val(2+(ii-1)*4:2+(ii-1)*4+2) = [0,0,0];
        p.domain_idx(2+(ii-1)*4:2+(ii-1)*4+2) = [ii,ii,ii];
        p.boundary_idx(2+(ii-1)*4:2+(ii-1)*4+2) = [1,1,1];
        p.element_idx(2+(ii-1)*4:2+(ii-1)*4+2) = [4,3,2]; 
    end
    
    p.xgeom(14:15) = [x0+6*tau_m,x0+6*tau_m];
    p.ygeom(14:15) = [y0,y0+h_c/2+g+h_m];
    p.domain_idx(14:15) = [4,11];
    p.boundary_idx(14:15) = [1,1];
    p.element_idx(14:15) = [40,13];
    bndr.idx(14:15) = [2,3];    

    for ii = 1:6
        p.xgeom(16+(ii-1)*4:16+(ii-1)*4+3) = [x0+(7-ii)*tau_m-dwm,x0+(7-ii)*tau_m-dwm, ...
            x0+(7-ii)*tau_m-dwm-w_m,x0+(7-ii)*tau_m-dwm-w_m];
        p.ygeom(16+(ii-1)*4:16+(ii-1)*4+3) = [y0+h_c/2+g+h_m,y0+h_c/2+g,y0+h_c/2+g,y0+h_c/2+g+h_m];     
        bndr.idx(16+(ii-1)*4:16+(ii-1)*4+3) = [3,3,3,3];
        bndr.val(16+(ii-1)*4:16+(ii-1)*4+3) = [0,mod(ii,2)*par.rho,0,0];
        p.domain_idx(16+(ii-1)*4:16+(ii-1)*4+3) = [4+ii,4+ii,4+ii,11];
        p.boundary_idx(16+(ii-1)*4:16+(ii-1)*4+3) = [1,1,1,1];
        p.element_idx(16+(ii-1)*4:16+(ii-1)*4+3) = [2,1,4,11-(ii-1)*2]; 
    end

    p.xgeom(40:41) = [x0,x0];
    p.ygeom(40:41) = [y0+h_c/2+g+h_m,y0];
    p.domain_idx(40) = [4];
    p.boundary_idx(40) = [1];
    p.element_idx(40) = [14];
    bndr.idx(40) = [2];    

    [Geom(4).boundary(1)] = createRandomBoundary(p,bndr,options);      

    %============================  Magnets  ===============================    
    for nn = 1:6
        p.x0 = x0+(7-nn)*tau_m-dwm-w_m;
        p.x1 = x0+(7-nn)*tau_m-dwm;
        p.y0 = y0+h_c/2+g;
        p.y1 = y0+h_c/2+g+h_m;
        bndr.idx = [3,3,3,3]; % idx = 1 is u bndr, idx = 0 is un bndr
        bndr.val = [mod(nn,2)*par.rho,0,0,0];
        p.domain_idx = [4,4,11,4]; 
        p.boundary_idx = [1,1,1,1];
        p.element_idx = [17+(nn-1)*4,16+(nn-1)*4,12-(nn-1)*2,18+(nn-1)*4];
        [Geom(4+nn).boundary(1)] = createSquareBoundary(p,bndr,options,0);          
    end
      
    %===========================  Back Iron  ==============================
    p.x0 = x0;
    p.y0 = y0;  

    p.xgeom = zeros(1,4+2*6+1);
    p.ygeom = zeros(1,4+2*6+1);
    p.domain_idx = zeros(1,4+2*6);
    p.boundary_idx = zeros(1,4+2*6);
    p.element_idx = zeros(1,4+2*6);
    bndr.idx = zeros(1,4+2*6);
    bndr.val = zeros(1,4+2*6);

    p.xgeom(1) = x0;
    p.ygeom(1) = y0+h_c/2+g+h_m;

    for nn = 1:6
        p.xgeom(2+(nn-1)*2:2+(nn-1)*2+1) = [x0+(nn-1)*tau_m+dwm,x0+(nn-1)*tau_m+dwm+w_m];
        p.ygeom(2+(nn-1)*2:2+(nn-1)*2+1) = [y0+h_c/2+g+h_m,y0+h_c/2+g+h_m];
        bndr.idx(1+(nn-1)*2:1+(nn-1)*2+1) = [3,3];
        bndr.val(1+(nn-1)*2:1+(nn-1)*2+1) = [0,0];
        p.domain_idx(1+(nn-1)*2:1+(nn-1)*2+1) = [4,4+7-nn];
        p.boundary_idx(1+(nn-1)*2:1+(nn-1)*2+1) = [1,1];
        p.element_idx(1+(nn-1)*2:1+(nn-1)*2+1) = [39-(nn-1)*4,3]; 
    end     
        p.xgeom(14:17) = [x0+6*tau_m,x0+6*tau_m,x0,x0];
        p.ygeom(14:17) = [y0+h_c/2+g+h_m,y0+h_c/2+g+h_m+h_bi,y0+h_c/2+g+h_m+h_bi,...
            y0+h_c/2+g+h_m];
        bndr.idx(13:16) = [3,2,1,2];    
        bndr.val(13:16) = [0,0,0,0];
        p.domain_idx(13:16) = [4,11,0,11];
        p.boundary_idx(13:16) = [1,1,0,1];
        p.element_idx(13:16) = [15,16,0,14]; 
        [Geom(11).boundary(1)] = createRandomBoundary(p,bndr,options);          
end