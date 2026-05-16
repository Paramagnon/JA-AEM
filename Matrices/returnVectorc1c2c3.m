function [c1,c2,c3] = returnVectorc1c2c3(Geom)
    % Find total matrix dimensions, all boundary elements (N) results in
    % NxN matrices for C1 and C2, and Nx1 C3 vector
    N = 0;
    for kk = 1:length(Geom)
        for ii = 1:length(Geom(kk).boundary)
            N = N + length(Geom(kk).boundary(ii).bndr_idx);
        end
    end
    c1 = zeros(N,1);       
    c2 = zeros(N,1);
    c3 = zeros(N,1);
    
    Noffset = 0;
    for kk = 1:length(Geom)
        for ll = 1:length(Geom(kk).boundary)
            Npart = length(Geom(kk).boundary(ll).c1);
            c1(1+Noffset:Noffset+Npart) = Geom(kk).boundary(ll).c1;
            c2(1+Noffset:Noffset+Npart) = Geom(kk).boundary(ll).c2;
            c3(1+Noffset:Noffset+Npart) = Geom(kk).boundary(ll).c3;
            Noffset = Noffset + Npart;
        end
    end
end