function [meshPointsX, meshPointsY, bndr_idx, bndr_idxT, Nidx, c1, c2, c3, c1T, c2T, c3T] = returnMeshCircleAEM(p,bndr,options,cw)
    % CW = 1 is clockwise, CW = 0 is anti-clockwise mesh definition
    x0 = p.x0;
    y0 = p.y0;
    r = p.r;
    dist = 2*pi*r;
    Nmin = 8;
    Nmesh = max(round(dist*options.N,0),Nmin);

    bndr_idx = ones(1,Nmesh-1)*bndr.idx;       
    bndr_idxT = ones(1,Nmesh-1)*bndr.idxT;       
    c1 = ones(Nmesh-1,1)*bndr.C1;
    c2 = ones(Nmesh-1,1)*bndr.C2;
    c3 = ones(Nmesh-1,1)*bndr.C3;
    c1T = ones(Nmesh-1,1)*bndr.C1T;
    c2T = ones(Nmesh-1,1)*bndr.C2T;
    c3T = ones(Nmesh-1,1)*bndr.C3T;

    theta = 0:2*pi/(Nmesh-1):2*pi;
    if cw
        theta = theta(end:-1:1);
    end
    
    meshPointsX = r*cos(theta) + x0;
    meshPointsY = r*sin(theta) + y0;
    Nidx = Nmesh-1;
end