function [C1T,C2T] = returnMatrixC1TC2T(Geom,GTlow)
    % Find total matrix dimensions, all boundary elements (N) results in
    % NxN matrices for C1 and C2, and Nx1 C3 vector
    
    N = 0;
    for kk = 1:length(Geom)
        for ii = 1:length(Geom(kk).boundary)
            N = N + length(Geom(kk).boundary(ii).bndr_idxT);
        end
    end
    C1T = zeros(N,N);
    C2T = zeros(N,N);    
    % Each domain will have an alpha, beta, and gamma per edge
    % This will be made into an C1 and C2 matrix the size of the boundary
    % of the domain in the meshing program
    
    % Follows the logic c1*u + c2*q = c3
    % Stored in geom as c1 = alpha, c2 = beta, c3 = gamma
    % Geom stores per edge an alpha, beta, gamma, and during meshing this
    % is created in to a c1 and c2 matrix and c3 vector per geometry
    
    % Location of boundary conditions in the matrix stay the same, the
    % values might change, so it is setup such that you get an index matrix
    % and a C1, C2 vector which remains the same. The updated values can
    % then directly be calculated
    
    Noffset = 0;
    for kk = 1:length(Geom)
        Npart = Geom(kk).Nbnd;
        C1geom = zeros(Npart,N);
        C2geom = zeros(Npart,N);       

        for ii = 1:length(Geom(kk).bndr_idxT)
            C1row = zeros(1,N);
            C2row = zeros(1,N);
            
            if Geom(kk).bndr_idxT(ii) >= 2
                Geomidx = Geom(kk).idx_linked(2,ii)-GTlow+1;
                Geomoffset = sum([Geom(1:Geomidx-1).Nbnd]);
                if Geomidx > kk % kk is C1 u, Geomidx is C2 q
                    C1row(Noffset+ii) = 1;
                    C1row(Geomoffset+Geom(kk).idx_linked(1,ii)) = 1;                    
                elseif Geomidx == kk
                    if Geom(kk).idx_linked(1,ii) > ii % (kk,ii) is C1 u, Geomidx is C2 q
                        C1row(Noffset+ii) = 1;
                        C1row(Noffset+Geom(kk).idx_linked(1,ii)) = 1;                        
                    else % linked is smaller than ii, (kk,ii) is C2 q, Geomidx is C1 u
                        C2row(Noffset+ii) = 1;
                        C2row(Noffset+Geom(kk).idx_linked(1,ii)) = 1;
                    end
                else
                    C2row(Noffset+ii) = 1;
                    C2row(Geomoffset+Geom(kk).idx_linked(1,ii)) = 1;                    
                end
            else
                C1row(Noffset+ii) = 1;
                C2row(Noffset+ii) = 1;
            end

            C1geom(ii,:) = C1row;
            C2geom(ii,:) = C2row;            
        end

        C1T(1+Noffset:Npart+Noffset,:) = C1geom;
        C2T(1+Noffset:Npart+Noffset,:) = C2geom;        
        Noffset = Noffset + Npart;
    end
    C1T = sparse(C1T);
    C2T = sparse(C2T);
end