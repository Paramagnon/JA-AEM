function [G,H] = calculateMatrixGHMultiBoundary(Geom)
    N = 0;   
    Nb = length(Geom.boundary);
    Nii = zeros(1,Nb+1);

    for ii = 1:Nb
        Nbnd = length(Geom.boundary(ii).Xm);
        Nii(ii+1) = Nii(ii) + Nbnd - 1; 
        N = N + Nbnd - 1;
    end

    G = zeros(N,N);
    H = zeros(N,N);
    
    for ii = 1:Nb
        X = Geom.boundary(ii).Xm;
        Y = Geom.boundary(ii).Ym;
        
        [G(Nii(ii)+1:Nii(ii+1),Nii(ii)+1:Nii(ii+1)),H(Nii(ii)+1:Nii(ii+1),Nii(ii)+1:Nii(ii+1))] = returnMatrixGHMultiDiagonal(X,Y);
    end
    if Nb > 1
        for nn = 1:Nb-1
            Nidx = nn;
            [G1,H1,G2,H2] = returnMatrixGHMultiOffDiagonal(Geom,Nidx,Nii);
            
            G(Nii(nn)+1:Nii(nn+1),Nii(nn+1)+1:end) = G1;
            G(Nii(nn+1)+1:end,Nii(nn)+1:Nii(nn+1)) = G2;
    
            H(Nii(nn)+1:Nii(nn+1),Nii(nn+1)+1:end) = H1;
            H(Nii(nn+1)+1:end,Nii(nn)+1:Nii(nn+1)) = H2;
        end
    end
end