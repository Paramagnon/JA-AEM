function [T,Qi,Ui] = returnInternalUiQiTThermal(Geom, SimGeom)
    %-----------------------------------------------
    % Approximation function used: phi = 1 + r 
    % Corresponding particular solution: u = 1/4 r^2 + 1/9 r^3

    % q = d u / d n (with respect to the boundary elements N)
    % = (du/ dr)*(dr/dn) = (1/2 r + 1/3 r^2)
    % r = sqrt( (ksi - x)^2 + (eta - y)^2 ) where ksi and eta are the
    % locations on the boundary and x and y the location of the nodes L
    % dr/dn = cos(beta)*(ksi-x)/r + sin(beta)*(eta-y)/r
    % where cos(beta) = d ksi/sqrt(d ksi^2 + d eta^2) and 
    % sin(beta) = d eta/sqrt(d ksi^2 + d eta^2) representing the normal
    % pointing away from the boundary element

    Trng = SimGeom.Rngs.Nrng;
    % ---------------------------------------------------------------------
    % Qi and Ui section
    % ---------------------------------------------------------------------          
    
    % Find total matrix dimensions    
    N = 0;
    L = 0;
    for kk = 1:length(Geom)
        N = N + Geom(kk).Nbnd;
        L = L + Geom(kk).LdmnT;
    end

    Ui = zeros(N,L);       
    Qi = zeros(N,L);
    
    GeomOffsetN = 0;
    GeomOffsetL = 0;
    for pp = 1:length(Geom)        
        Npart = Geom(pp).Nbnd;      
        Lpart = Geom(pp).LdmnT;

        Uip = zeros(Npart,Lpart);
        Qip = zeros(Npart,Lpart);
       
        epsilon = SimGeom.dmnrT(pp);
        for ii = 1:Geom(pp).LdmnT
            x = Geom(pp).domainT.Xn(ii); % Internal node (y)
            y = Geom(pp).domainT.Yn(ii); % Internal node (y)
            
            Noffset = 0;
            for kk = 1:length(Geom(pp).boundary)
                Xn = Geom(pp).boundary(kk).Xn;
                Yn = Geom(pp).boundary(kk).Yn;
    
                Xm = Geom(pp).boundary(kk).Xm;
                Ym = Geom(pp).boundary(kk).Ym;
    
                for mm = 1:length(Geom(pp).boundary(kk).Xn)
                    ksi = Xn(mm); % Boundary nodes (x)
                    eta = Yn(mm); % Boundary nodes (y)
    
                    r = sqrt((ksi-x)^2+(eta-y)^2);
                    if r < 1e-12
                        Uip(Noffset+mm,ii) = 0;
                        Qip(Noffset+mm,ii) = 0;
                    else                       
                        reps = r/epsilon;
                        r2 = reps^2;
                        Uip(Noffset+mm,ii) = 0.25*r2 + (1/9)*r2*reps;
                        dksi = Xm(mm+1) - Xm(mm);
                        deta = Ym(mm+1) - Ym(mm);
                        sqrdr = 1/(sqrt(dksi^2+deta^2)*r);
                        % nx = dym/sqrt(dxm^2+dym^2);
                        % ny = -dxm/sqrt(dxm^2+dym^2);
                        drdn = ((ksi-x)*deta - (eta-y)*dksi)*sqrdr;
                        Qip(Noffset+mm,ii) = (0.5*reps + (1/3)*r2)*drdn;                                               
                    end
                end
                Noffset = Noffset + length(Geom(pp).boundary(kk).Xn);
            end
        end
        Ui(1+GeomOffsetN:Npart+GeomOffsetN,1+GeomOffsetL:Lpart+GeomOffsetL) = Uip;
        Qi(1+GeomOffsetN:Npart+GeomOffsetN,1+GeomOffsetL:Lpart+GeomOffsetL) = Qip;

        GeomOffsetN = GeomOffsetN + Npart;
        GeomOffsetL = GeomOffsetL + Lpart;
    end
 
    % ---------------------------------------------------------------------
    % T matrix: T = (H Ui - G Qi)
    % ---------------------------------------------------------------------

    if nnz(Ui)/numel(Ui) < 0.1   
        Ui = sparse(Ui);
        Qi = sparse(Qi);
    end  
    T = (SimGeom.H(Trng,Trng) * Ui - SimGeom.G(Trng,Trng) * Qi);
end