function [mu] = returnPropertiesMagnetics(Geom,SimGeom)

    Ltot = sum([Geom.Ldmn]);
    % TT = SimGeom.TLT;
    mu = zeros(Ltot,1);

    Loffset = 0; 
    Toffset = 0;   
    for ii = 1:length(Geom)
        Lpart = Geom(ii).Ldmn;
        Tpart = Geom(ii).LdmnT;       
        
        if Lpart == Tpart
            mu(Loffset+1:Loffset+Lpart) = Geom(ii).mu_f(TT(Toffset+1:Toffset+Tpart));
        end

        Loffset = Loffset + Lpart;
        Toffset = Toffset + Tpart;        
    end   
end