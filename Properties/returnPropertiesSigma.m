function [sigmaA] = returnPropertiesSigma(Geom,SimGeom)

    Ltot = sum([Geom.Ldmn]);
    % Ttot = sum([Geom.LdmnT]);
    % Mtot = sum([SimGeom.PlotGeom.Ldmn]);
    
    TT = SimGeom.TLT;
    % TPlt = SimGeom.TLPlt;

    sigmaA = zeros(Ltot,1);
    % sigmaT = zeros(Ttot,1);
    % sigmaPlt = zeros(Mtot,1);

    Loffset = 0; 
    Toffset = 0;
    % Moffset = 0;
    for ii = 1:length(Geom)
        Lpart = Geom(ii).Ldmn;
        Tpart = Geom(ii).LdmnT;
        % Mpart = SimGeom.PlotGeom(ii).Ldmn;
        
        % sigmaT(Toffset+1:Toffset+Tpart) = Geom(ii).sigma_f(TT(Toffset+1:Toffset+Tpart));
        % sigmaPlt(Moffset+1:Moffset+Mpart) = Geom(ii).sigma_f(TPlt(Moffset+1:Moffset+Mpart));

        if Lpart == Tpart
            % sigmaA(Loffset+1:Loffset+Lpart) = sigmaT(Toffset+1:Toffset+Tpart);
            sigmaA(Loffset+1:Loffset+Lpart) = Geom(ii).sigma_f(TT(Toffset+1:Toffset+Tpart));           
        end

        Loffset = Loffset + Lpart;
        Toffset = Toffset + Tpart;
        % Moffset = Moffset + Mpart;
    end   

end