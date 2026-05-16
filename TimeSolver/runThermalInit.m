function ThermalGeom = runThermalInit(Geom,par,options)
    TL = sum(Geom(1:options.TgeomL-1).Nbnd) + 1;
    TU = sum(Geom(1:options.TgeomU.Nbnd));

    Geom(TL:TU)
    
    for kk = 1:length(Geom)
        [Geom(kk).G,Geom(kk).H] = calculateMatrixGHMultiBoundary(Geom(kk));
        [Geom(kk).bndr_idx] = returnMultiDomainBndrIdx(Geom(kk));
    end

    [ThermalGeom.G00plt,ThermalGeom.H00plt,ThermalGeom.U00plt,ThermalGeom.PlotGeom] = returnExportMatrices(Geom, options);
    [ThermalGeom.C1idx,ThermalGeom.C2idx] = returnMatrixC1C2(Geom);
    [ThermalGeom.G,ThermalGeom.H] = returnMatrixGH(Geom);
    [ThermalGeom.T,ThermalGeom.Qi,ThermalGeom.Ui] = returnInternalUiQiT(Geom,ThermalGeom);
    [ThermalGeom.G00,ThermalGeom.Gx0,ThermalGeom.Gy0,ThermalGeom.H00,ThermalGeom.Hx0,ThermalGeom.Hy0,ThermalGeom.U00,ThermalGeom.Ux0,ThermalGeom.Uy0,ThermalGeom.Uxx,ThermalGeom.Uyy] = returnAlphaMatrices(Geom);
    ThermalGeom.Temp = par.T0*ones(sum([Geom.Ldmn]),1);
    ThermalGeom.TempPlt = par.T0*ones(sum([ThermalGeom.PlotGeom.Ldmn]),1);
    % Split function into mu and sigma
    [ThermalGeom.sigma, ThermalGeom.mu] = returnPropertiesMagnetics(Geom,ThermalGeom);
    [ThermalGeom.sigma, ThermalGeom.sigmaplt] = returnPropertiesSigma(Geom,ThermalGeom);
    [ThermalGeom.an, ThermalGeom.bn] = boundaryFourierInnit(par);
end