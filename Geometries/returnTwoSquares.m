function Geom = returnTwoSquares(options)
    Geom(1).a = 1; % Domain property
    Geom(2).a = 1; % Domain property
    Geom(3).a = 1; % Domain property
    Geom(1).source = 0;
    Geom(2).source = 0;
    Geom(3).source = 0;

    % --- Outer square
    Geom(1).boundary(1).x0 = 0;
    Geom(1).boundary(1).x1 = 5;
    Geom(1).boundary(1).y0 = 0;
    Geom(1).boundary(1).y1 = 5;
    % --- Inner square 1 
    Geom(2).boundary(1).x0 = 5;
    Geom(2).boundary(1).x1 = 10;
    Geom(2).boundary(1).y0 = 0;
    Geom(2).boundary(1).y1 = 5;

    % --- Inner square 3 test
    Geom(3).boundary(1).x0 = 0;
    Geom(3).boundary(1).x1 = 5;
    Geom(3).boundary(1).y0 = 5;
    Geom(3).boundary(1).y1 = 10;


    bndr.idx = [1,3,3,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,0,0,40];
    Geom(1).boundary(1).domain_idx = [0,2,3,0]; %[0,2,3,0];
    Geom(1).boundary(1).boundary_idx = [0,1,1,0];
    Geom(1).boundary(1).element_idx = [0,4,1,0];
    [Geom(1).boundary(1).xgeom, Geom(1).boundary(1).ygeom] = returnSquarePoints(Geom(1).boundary(1),0);            
    [Geom(1).boundary(1).Xm, Geom(1).boundary(1).Ym, Geom(1).boundary(1).bndr_idx, Geom(1).boundary(1).bndr_val, Geom(1).boundary(1).Nidx] = returnGeomPoints(Geom(1).boundary(1),bndr,options);
    

    bndr.idx = [1,0,1,3]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,-5,40,0];
    Geom(2).boundary(1).domain_idx = [0,0,0,1];
    Geom(2).boundary(1).boundary_idx = [0,0,0,1];
    Geom(2).boundary(1).element_idx = [0,0,0,2];   
    [Geom(2).boundary(1).xgeom, Geom(2).boundary(1).ygeom] = returnSquarePoints(Geom(2).boundary(1),0);        
    [Geom(2).boundary(1).Xm, Geom(2).boundary(1).Ym, Geom(2).boundary(1).bndr_idx, Geom(2).boundary(1).bndr_val, Geom(2).boundary(1).Nidx] = returnGeomPoints(Geom(2).boundary(1),bndr,options);


    bndr.idx = [3,0,0,1]; % idx = 1 is u bndr, idx = 0 is un bndr
    bndr.val = [0,-10,-10,40];
    Geom(3).boundary(1).domain_idx = [1,0,0,0];
    Geom(3).boundary(1).boundary_idx = [1,0,0,0];
    Geom(3).boundary(1).element_idx = [3,0,0,0];   
    [Geom(3).boundary(1).xgeom, Geom(3).boundary(1).ygeom] = returnSquarePoints(Geom(3).boundary(1),0);        
    [Geom(3).boundary(1).Xm, Geom(3).boundary(1).Ym, Geom(3).boundary(1).bndr_idx, Geom(3).boundary(1).bndr_val, Geom(3).boundary(1).Nidx] = returnGeomPoints(Geom(3).boundary(1),bndr,options);
end