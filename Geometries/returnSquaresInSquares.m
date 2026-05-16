function Geom = returnSquaresInSquares(options)
    Geom(1).a = 0;
    Geom(1).source = 0;
    % --- Outer square
    Geom(1).boundary(1).x0 = 0;
    Geom(1).boundary(1).x1 = 5;
    Geom(1).boundary(1).y0 = 0;
    Geom(1).boundary(1).y1 = 5;
    % --- Inner square 1 
    Geom(1).boundary(2).x0 = 3.75;
    Geom(1).boundary(2).x1 = 4.25;
    Geom(1).boundary(2).y0 = 3.5;
    Geom(1).boundary(2).y1 = 4;
    % --- Inner square 2
    Geom(1).boundary(3).x0 = 0.75;
    Geom(1).boundary(3).x1 = 1.25;
    Geom(1).boundary(3).y0 = 3.5;
    Geom(1).boundary(3).y1 = 4;
    % --- Inner square 3
    Geom(1).boundary(4).x0 = 2.25;
    Geom(1).boundary(4).x1 = 2.75;
    Geom(1).boundary(4).y0 = 0.75;
    Geom(1).boundary(4).y1 = 1.25;
    
    % Geom(1).boundary(6).x0 = 0.25;
    % Geom(1).boundary(6).x1 = 1;
    % Geom(1).boundary(6).y0 = 0.25;
    % Geom(1).boundary(6).y1 = 1;
        
    bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [10,10,10,10];
    [Geom(1).boundary(1).xgeom, Geom(1).boundary(1).ygeom] = returnSquarePoints(Geom(1).boundary(1),0);
    [Geom(1).boundary(1).Xm, Geom(1).boundary(1).Ym, Geom(1).boundary(1).bndr_idx, Geom(1).boundary(1).bndr_val, Geom(1).boundary(1).Nidx] = returnGeomPoints(Geom(1).boundary(1),bndr,options);
    
    bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [-50,-50,-50,-50];
    [Geom(1).boundary(2).xgeom, Geom(1).boundary(2).ygeom] = returnSquarePoints(Geom(1).boundary(2),1);       
    [Geom(1).boundary(2).Xm, Geom(1).boundary(2).Ym, Geom(1).boundary(2).bndr_idx, Geom(1).boundary(2).bndr_val, Geom(1).boundary(2).Nidx] = returnGeomPoints(Geom(1).boundary(2),bndr,options);
    
    bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [-50,-50,-50,-50];
    [Geom(1).boundary(3).xgeom, Geom(1).boundary(3).ygeom] = returnSquarePoints(Geom(1).boundary(3),1);        
    [Geom(1).boundary(3).Xm, Geom(1).boundary(3).Ym, Geom(1).boundary(3).bndr_idx, Geom(1).boundary(3).bndr_val, Geom(1).boundary(3).Nidx] = returnGeomPoints(Geom(1).boundary(3),bndr,options);
    
    bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [-50,-50,-50,-50];
    [Geom(1).boundary(4).xgeom, Geom(1).boundary(4).ygeom] = returnSquarePoints(Geom(1).boundary(4),1);       
    [Geom(1).boundary(4).Xm, Geom(1).boundary(4).Ym, Geom(1).boundary(4).bndr_idx, Geom(1).boundary(4).bndr_val, Geom(1).boundary(4).Nidx] = returnGeomPoints(Geom(1).boundary(4),bndr,options);
      
    Geom(1).boundary(5).p.x0 = 2.5;
    Geom(1).boundary(5).p.y0 = 2.5;
    Geom(1).boundary(5).p.r = 0.5;
    N = 60;
    bndr.idx = [3];
    bndr.val = [100];
    [Geom(1).boundary(5).Xm, Geom(1).boundary(5).Ym, Geom(1).boundary(5).bndr_idx, Geom(1).boundary(5).bndr_val, Geom(1).boundary(5).Nidx] = returnMeshCircle(Geom(1).boundary(5).p,bndr,N,1);
    Geom(1).boundary(5).xgeom = [Geom(1).boundary(5).Xm(1) Geom(1).boundary(5).Xm(2)];
    Geom(1).boundary(5).ygeom = [Geom(1).boundary(5).Ym(1) Geom(1).boundary(5).Ym(2)];

    Geom(2).boundary(1).p.x0 = 2.5;
    Geom(2).boundary(1).p.y0 = 2.5;
    Geom(2).boundary(1).p.r = 0.5;
    N = 60;
    bndr.idx = [1];
    bndr.val = [100];
    [Geom(2).boundary(1).Xm, Geom(2).boundary(1).Ym, Geom(2).boundary(1).bndr_idx, Geom(2).boundary(1).bndr_val, Geom(2).boundary(1).Nidx] = returnMeshCircle(Geom(2).boundary(1).p,bndr,N,0);
    Geom(2).boundary(1).xgeom = [Geom(2).boundary(1).Xm(1) Geom(2).boundary(1).Xm(2)];
    Geom(2).boundary(1).ygeom = [Geom(2).boundary(1).Ym(1) Geom(2).boundary(1).Ym(2)];

    % bndr.idx = [1,1,1,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    % bndr.val = [50,50,50,50];
    % [Geom(1).boundary(6).xgeom, Geom(1).boundary(6).ygeom] = returnSquarePoints(Geom(1).boundary(6),1);       
    % [Geom(1).boundary(6).Xm, Geom(1).boundary(6).Ym, Geom(1).boundary(6).bndr_idx, Geom(1).boundary(6).bndr_val, Geom(1).boundary(6).Nidx] = returnGeomPoints(Geom(1).boundary(6),bndr,options);
    Geom(1).domain = insertDomainPoints(Geom, options);
end