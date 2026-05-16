function [bndr_idx,bndr_idxT] = returnMultiDomainBndrIdx(Geom)
    N = 0;   
    Nb = length(Geom.boundary);
    for ii = 1:Nb
        Nbnd = length(Geom.boundary(ii).Xm);        
        N = N + Nbnd - 1;
    end
    bndr_idx = zeros(1,N);
    bndr_idxT = zeros(1,N);

    Noffset = 0;
    for ii = 1:Nb
        Nbnd = length(Geom.boundary(ii).Xn);        
        bndr_idx(1+Noffset:Nbnd+Noffset) = Geom.boundary(ii).bndr_idx;
        % bndr_idxT(1+Noffset:Nbnd+Noffset) = Geom.boundary(ii).bndr_idxT;
        Noffset = Noffset+Nbnd;
    end

end