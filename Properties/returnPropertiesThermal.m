function [kappa, rho, Cp] = returnPropertiesThermal(Geom,SimGeom)
    Ttot = sum([Geom.LdmnT]);
    TT = SimGeom.TLT;

    Cp = zeros(Ttot,1);
    kappa = zeros(Ttot,1);
    rho = zeros(Ttot,1);   

    Toffset = 0;
    for ii = 1:length(Geom)
        Tpart = Geom(ii).LdmnT;
        kappa(Toffset+1:Toffset+Tpart) = Geom(ii).kappa_f(TT(Toffset+1:Toffset+Tpart),Geom(ii).kappa_typ);
        rho(Toffset+1:Toffset+Tpart) = Geom(ii).rho_f(TT(Toffset+1:Toffset+Tpart));
        Cp(Toffset+1:Toffset+Tpart) = Geom(ii).Cp_f(TT(Toffset+1:Toffset+Tpart));
        Toffset = Toffset + Tpart;       
    end    
end