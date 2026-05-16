function [meshPointsX, meshPointsY, bndr_idx, bndr_idxT, Nidx, C1, C2, C3, C1T, C2T, C3T] = returnGeomPointsAEM(p,bndr,options)
    % Define p.x0, p.y0, p.xgeom, p.ygeom 
    % Returns meshpoints and boundary values and type
    dmin = options.dmin;
    gamma = options.gamma;
    N = options.N;
    Nx = options.Nx;
    Ny = options.Ny;
    Nmin = options.Nmin;
    linear = options.linear; 
    
    xgeom = p.xgeom; 
    ygeom = p.ygeom;
      
    meshPointsX = [xgeom(1)];
    meshPointsY = [ygeom(1)];
    
    bndr_idx = [];    
    bndr_idxT = [];    

    C1 = [];
    C2 = [];
    C3 = [];

    C1T = [];
    C2T = [];
    C3T = [];

    Nidx = zeros(1,length(xgeom)-1);
    if linear  
        dx = (xgeom(1:end-1)-xgeom(2:end));
        dy = (ygeom(1:end-1)-ygeom(2:end));        
        dist = sqrt(dx.^2+dy.^2);

        Nsplit = zeros(size(dist));
        for kk = 1:length(xgeom)-1
            if dx(kk) == 0 
                Nsplit(kk) = round(Ny*dist(kk),0);
            elseif dy(kk) == 0
                Nsplit(kk) = round(Nx*dist(kk),0);
            else
                Nsplit(kk) = round(N*dist(kk),0);
            end
            if Nsplit(kk) < Nmin
                Nsplit(kk) = Nmin;
            end
            % Nsplit = N*dist;
            % Nsplit = round(Nsplit,0);
            % N specifies mesh elements per length unit
        end

        for kk = 1:length(xgeom)-1
            [Xint, Yint] = returnMeshLinearAEM(xgeom(kk),ygeom(kk),xgeom(kk+1),ygeom(kk+1),Nsplit(kk));
            meshPointsX(end+1:end+length(Xint)-1) = Xint(2:end);
            meshPointsY(end+1:end+length(Yint)-1) = Yint(2:end);    
            
            bndr_idx(end+1:end+length(Xint)-1) = bndr.idx(kk);    
            bndr_idxT(end+1:end+length(Xint)-1) = bndr.idxT(kk);    

            C1(end+1:end+length(Xint)-1) = bndr.C1(kk);
            C2(end+1:end+length(Xint)-1) = bndr.C2(kk);
            C3(end+1:end+length(Xint)-1) = bndr.C3(kk);
            C1T(end+1:end+length(Xint)-1) = bndr.C1T(kk);
            C2T(end+1:end+length(Xint)-1) = bndr.C2T(kk);
            C3T(end+1:end+length(Xint)-1) = bndr.C3T(kk);
            
            Nidx(kk) = length(Xint)-1;
        end
        bndr_idx = bndr_idx(1:end);        
        bndr_idxT = bndr_idxT(1:end);        
    else
        
        for kk = 1:length(xgeom)-1
            [Xint, Yint] = returnMeshExponentialAEM(xgeom(kk),ygeom(kk),xgeom(kk+1),ygeom(kk+1),dmin,gamma);
            meshPointsX(end+1:end+length(Xint)-1) = Xint(2:end);
            meshPointsY(end+1:end+length(Yint)-1) = Yint(2:end);    
            
            bndr_idx(end+1:end+length(Xint)-1) = bndr.idx(kk);            
            bndr_idxT(end+1:end+length(Xint)-1) = bndr.idxT(kk);            

            C1(end+1:end+length(Xint)-1) = bndr.C1(kk);
            C2(end+1:end+length(Xint)-1) = bndr.C2(kk);
            C3(end+1:end+length(Xint)-1) = bndr.C3(kk);
            C1T(end+1:end+length(Xint)-1) = bndr.C1T(kk);
            C2T(end+1:end+length(Xint)-1) = bndr.C2T(kk);
            C3T(end+1:end+length(Xint)-1) = bndr.C3T(kk);

            Nidx(kk) = length(Xint)-1;
        end
        bndr_idx = bndr_idx(1:end);
        bndr_idxT = bndr_idxT(1:end);
    end
end
