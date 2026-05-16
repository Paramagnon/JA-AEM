function [G,H] = returnMatrixGH(Geom)
    % Simulation is setup such that the shared boundaries are handled in
    % the C1 and C2 matrix, so here it is just a simple insertion of each
    % square matrix for each Geom element on the diagonal
    % Each Geom satisfies the boundary equation H u - G q = T alpha

    % Find total matrix dimensions    
    N = 0;
    for kk = 1:length(Geom)
        for ii = 1:length(Geom(kk).boundary)
            N = N + length(Geom(kk).boundary(ii).bndr_idx);
        end
    end

    G = zeros(N,N);       
    H = zeros(N,N);
    
    Noffset = 0;
    for kk = 1:length(Geom)
        Npart = Geom(kk).Nbnd;      
        
        G(1+Noffset:Npart+Noffset,1+Noffset:Npart+Noffset) = Geom(kk).G;
        H(1+Noffset:Npart+Noffset,1+Noffset:Npart+Noffset) = Geom(kk).H;

        Noffset = Noffset + Npart;
    end  
    if nnz(G)/numel(G) < 0.1   
        G = sparse(G);
        H = sparse(H);
    end   
end