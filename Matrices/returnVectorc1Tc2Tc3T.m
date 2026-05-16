function [c1T,c2T,c3T] = returnVectorc1Tc2Tc3T(Geom)
    % Find total matrix dimensions, all boundary elements (N) results in
    % NxN matrices for C1 and C2, and Nx1 C3 vector
    N = 0;
    for kk = 1:length(Geom)
        for ii = 1:length(Geom(kk).boundary)
            N = N + length(Geom(kk).boundary(ii).bndr_idx);
        end
    end
    c1T = zeros(N,1);       
    c2T = zeros(N,1);
    c3T = zeros(N,1);
    
    Noffset = 0;
    for kk = 1:length(Geom)
        for ll = 1:length(Geom(kk).boundary)
            Npart = length(Geom(kk).boundary(ll).c1T);
            c1T(1+Noffset:Noffset+Npart) = Geom(kk).boundary(ll).c1T;
            c2T(1+Noffset:Noffset+Npart) = Geom(kk).boundary(ll).c2T;
            c3T(1+Noffset:Noffset+Npart) = Geom(kk).boundary(ll).c3T;
            Noffset = Noffset + Npart;
        end
    end
end