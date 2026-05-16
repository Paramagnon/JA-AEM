function [Geom] = returnThreePhaseCable(options,p)  
    N = 60;
    eps0 =  8.8541878188e-12;

    %--- Wire 1 - Vu, Iu
    Geom(1).a = 5*eps0;
    % Wire 1, inner circle
    Geom(1).boundary(1).p.x0 = p.x0 + p.L/2;
    Geom(1).boundary(1).p.y0 = p.y0 + p.L/2 + p.d_w;
    Geom(1).boundary(1).p.r = p.r_w;  
    Geom(1).boundary(1).domain_idx = [0];    
    bndr.idx = [1];
    bndr.val = [p.Vu];
    [Geom(1).boundary(1).Xm, Geom(1).boundary(1).Ym, Geom(1).boundary(1).bndr_idx, Geom(1).boundary(1).bndr_val, Geom(1).boundary(1).Nidx] = returnMeshCircle(Geom(1).boundary(1).p,bndr,N,1);
    Geom(1).boundary(1).xgeom = [Geom(1).boundary(1).Xm(1) Geom(1).boundary(1).Xm(2)];
    Geom(1).boundary(1).ygeom = [Geom(1).boundary(1).Ym(1) Geom(1).boundary(1).Ym(2)];
    % Wire 1, outer circle
    Geom(1).boundary(2).p.x0 = p.x0 + p.L/2;
    Geom(1).boundary(2).p.y0 = p.y0 + p.L/2 + p.d_w;
    Geom(1).boundary(2).p.r = p.r_w + p.d_ins; 
    Geom(1).boundary(2).domain_idx = [4];
    Geom(1).boundary(2).boundary_idx = [2];
    Geom(1).boundary(2).element_idx = [1];
    bndr.idx = [3];
    bndr.val = [0];
    [Geom(1).boundary(2).Xm, Geom(1).boundary(2).Ym, Geom(1).boundary(2).bndr_idx, Geom(1).boundary(2).bndr_val, Geom(1).boundary(2).Nidx] = returnMeshCircle(Geom(1).boundary(2).p,bndr,N,0);
    Geom(1).boundary(2).xgeom = [Geom(1).boundary(2).Xm(1) Geom(1).boundary(2).Xm(2)];
    Geom(1).boundary(2).ygeom = [Geom(1).boundary(2).Ym(1) Geom(1).boundary(2).Ym(2)];

    %--- Wire 2 - Vv, Iv
    Geom(2).a = 5*eps0;
    % Wire 2, inner circle
    Geom(2).boundary(1).p.x0 = p.x0 + p.L/2 + cos(-pi/6)*p.d_w;
    Geom(2).boundary(1).p.y0 = p.y0 + p.L/2 + sin(-pi/6)*p.d_w;
    Geom(2).boundary(1).p.r = p.r_w;  
    Geom(2).boundary(1).domain_idx = [0];    
    bndr.idx = [1];
    bndr.val = [p.Vv];
    [Geom(2).boundary(1).Xm, Geom(2).boundary(1).Ym, Geom(2).boundary(1).bndr_idx, Geom(2).boundary(1).bndr_val, Geom(2).boundary(1).Nidx] = returnMeshCircle(Geom(2).boundary(1).p,bndr,N,1);
    Geom(2).boundary(1).xgeom = [Geom(2).boundary(1).Xm(1) Geom(2).boundary(1).Xm(2)];
    Geom(2).boundary(1).ygeom = [Geom(2).boundary(1).Ym(1) Geom(2).boundary(1).Ym(2)];
    % Wire 2, outer circle
    Geom(2).boundary(2).p.x0 = p.x0 + p.L/2 + cos(-pi/6)*p.d_w;
    Geom(2).boundary(2).p.y0 = p.y0 + p.L/2 + sin(-pi/6)*p.d_w;
    Geom(2).boundary(2).p.r = p.r_w + p.d_ins;  
    Geom(2).boundary(2).domain_idx = [4];
    Geom(2).boundary(2).boundary_idx = [3];
    Geom(2).boundary(2).element_idx = [1];
    bndr.idx = [3];
    bndr.val = [0];
    [Geom(2).boundary(2).Xm, Geom(2).boundary(2).Ym, Geom(2).boundary(2).bndr_idx, Geom(2).boundary(2).bndr_val, Geom(2).boundary(2).Nidx] = returnMeshCircle(Geom(2).boundary(2).p,bndr,N,0);
    Geom(2).boundary(2).xgeom = [Geom(2).boundary(2).Xm(1) Geom(2).boundary(2).Xm(2)];
    Geom(2).boundary(2).ygeom = [Geom(2).boundary(2).Ym(1) Geom(2).boundary(2).Ym(2)];
    
    %--- Wire 2 - Vw, Iw
    Geom(3).a = 5*eps0;
    % Wire 3, inner circle
    Geom(3).boundary(1).p.x0 = p.x0 + p.L/2 + cos(-5*pi/6)*p.d_w;
    Geom(3).boundary(1).p.y0 = p.y0 + p.L/2 + sin(-5*pi/6)*p.d_w;
    Geom(3).boundary(1).p.r = p.r_w;  
    Geom(3).boundary(1).domain_idx = [0];
    bndr.idx = [1];
    bndr.val = [p.Vw];
    [Geom(3).boundary(1).Xm, Geom(3).boundary(1).Ym, Geom(3).boundary(1).bndr_idx, Geom(3).boundary(1).bndr_val, Geom(3).boundary(1).Nidx] = returnMeshCircle(Geom(3).boundary(1).p,bndr,N,1);
    Geom(3).boundary(1).xgeom = [Geom(3).boundary(1).Xm(1) Geom(3).boundary(1).Xm(2)];
    Geom(3).boundary(1).ygeom = [Geom(3).boundary(1).Ym(1) Geom(3).boundary(1).Ym(2)];
    % Wire 3, outer circle
    Geom(3).boundary(2).p.x0 = p.x0 + p.L/2 + cos(-5*pi/6)*p.d_w;
    Geom(3).boundary(2).p.y0 = p.y0 + p.L/2 + sin(-5*pi/6)*p.d_w;
    Geom(3).boundary(2).p.r = p.r_w + p.d_ins;  
    Geom(3).boundary(2).domain_idx = [4];
    Geom(3).boundary(2).boundary_idx = [4];
    Geom(3).boundary(2).element_idx = [1];
    bndr.idx = [3];
    bndr.val = [0];
    [Geom(3).boundary(2).Xm, Geom(3).boundary(2).Ym, Geom(3).boundary(2).bndr_idx, Geom(3).boundary(2).bndr_val, Geom(3).boundary(2).Nidx] = returnMeshCircle(Geom(3).boundary(2).p,bndr,N,0);
    Geom(3).boundary(2).xgeom = [Geom(3).boundary(2).Xm(1) Geom(3).boundary(2).Xm(2)];
    Geom(3).boundary(2).ygeom = [Geom(3).boundary(2).Ym(1) Geom(3).boundary(2).Ym(2)];
    
    %--- Outside region
    Geom(4).a = 1*eps0;
    % Outside box
    Geom(4).boundary(1).x0 = p.x0;
    Geom(4).boundary(1).x1 = p.x0+p.L;
    Geom(4).boundary(1).y0 = p.y0;
    Geom(4).boundary(1).y1 = p.y0+p.L;

    bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,0,0];
    Geom(4).boundary(1).domain_idx = [0,0,0,0];
    Geom(4).boundary(1).boundary_idx = [0,0,0,0];
    Geom(4).boundary(1).element_idx = [0,0,0,0];
    [Geom(4).boundary(1).xgeom, Geom(4).boundary(1).ygeom] = returnSquarePoints(Geom(4).boundary(1),0);
    [Geom(4).boundary(1).Xm, Geom(4).boundary(1).Ym, Geom(4).boundary(1).bndr_idx, Geom(4).boundary(1).bndr_val, Geom(4).boundary(1).Nidx] = returnGeomPoints(Geom(4).boundary(1),bndr,options);
   
    % Wire 1, outer circle
    Geom(4).boundary(2).p.x0 = p.x0 + p.L/2;
    Geom(4).boundary(2).p.y0 = p.y0 + p.L/2 + p.d_w;
    Geom(4).boundary(2).p.r = p.r_w + p.d_ins;  
    Geom(4).boundary(2).domain_idx = [1];
    Geom(4).boundary(2).boundary_idx = [1];
    Geom(4).boundary(2).element_idx = [1];
    bndr.idx = [3];
    bndr.val = [0];
    [Geom(4).boundary(2).Xm, Geom(4).boundary(2).Ym, Geom(4).boundary(2).bndr_idx, Geom(4).boundary(2).bndr_val, Geom(4).boundary(2).Nidx] = returnMeshCircle(Geom(4).boundary(2).p,bndr,N,1);
    Geom(4).boundary(2).xgeom = [Geom(4).boundary(2).Xm(1) Geom(4).boundary(2).Xm(2)];
    Geom(4).boundary(2).ygeom = [Geom(4).boundary(2).Ym(1) Geom(4).boundary(2).Ym(2)];

    % Wire 2, outer circle
    Geom(4).boundary(3).p.x0 = p.x0 + p.L/2 + cos(-pi/6)*p.d_w;
    Geom(4).boundary(3).p.y0 = p.y0 + p.L/2 + sin(-pi/6)*p.d_w;
    Geom(4).boundary(3).p.r = p.r_w + p.d_ins;  
    Geom(4).boundary(3).domain_idx = [2];
    Geom(4).boundary(3).boundary_idx = [1];
    Geom(4).boundary(3).element_idx = [1];
    bndr.idx = [3];
    bndr.val = [0];
    [Geom(4).boundary(3).Xm, Geom(4).boundary(3).Ym, Geom(4).boundary(3).bndr_idx, Geom(4).boundary(3).bndr_val, Geom(4).boundary(3).Nidx] = returnMeshCircle(Geom(4).boundary(3).p,bndr,N,1);
    Geom(4).boundary(3).xgeom = [Geom(4).boundary(3).Xm(1) Geom(4).boundary(3).Xm(2)];
    Geom(4).boundary(3).ygeom = [Geom(4).boundary(3).Ym(1) Geom(4).boundary(3).Ym(2)];

    % Wire 3, outer circle
    Geom(4).boundary(4).p.x0 = p.x0 + p.L/2 + cos(-5*pi/6)*p.d_w;
    Geom(4).boundary(4).p.y0 = p.y0 + p.L/2 + sin(-5*pi/6)*p.d_w;
    Geom(4).boundary(4).p.r = p.r_w + p.d_ins;  
    Geom(4).boundary(4).domain_idx = [3];
    Geom(4).boundary(4).boundary_idx = [1];
    Geom(4).boundary(4).element_idx = [1];
    bndr.idx = [3];
    bndr.val = [0];
    [Geom(4).boundary(4).Xm, Geom(4).boundary(4).Ym, Geom(4).boundary(4).bndr_idx, Geom(4).boundary(4).bndr_val, Geom(4).boundary(4).Nidx] = returnMeshCircle(Geom(4).boundary(4).p,bndr,N,1);
    Geom(4).boundary(4).xgeom = [Geom(4).boundary(4).Xm(1) Geom(4).boundary(4).Xm(2)];
    Geom(4).boundary(4).ygeom = [Geom(4).boundary(4).Ym(1) Geom(4).boundary(4).Ym(2)];
end

