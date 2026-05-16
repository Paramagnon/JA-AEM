function SimSol = returnSimSol(SimGeom,alpha,alphaT,options)
    % rngN = SimGeom.Rngs.Nrng;
    % rngM = SimGeom.Rngs.Mrng;

    %--- Magnetics
    S00 = [SimGeom.H00plt, -SimGeom.G00plt]*SimGeom.Z + SimGeom.U00plt + SimGeom.G00plt*SimGeom.Qi - SimGeom.H00plt*SimGeom.Ui;
    z00 = [SimGeom.H00plt, -SimGeom.G00plt]*SimGeom.e;   
    F00 = SimGeom.puL + SimGeom.G00plt*SimGeom.pqN - SimGeom.H00plt*SimGeom.puN;
    SimSol.uL = S00*alpha + z00 + F00;
    
    uq = SimGeom.Z*alpha + SimGeom.e;
    N = length(uq)/2;
    SimSol.uN = uq(1:N);
    SimSol.uknwn = [SimSol.uN; SimSol.uL];

    % Implement force per depth, use coil geom, make a Gx0 Gy0 matrice for
    % it and calculate force by contour integral


    %--- Thermal
    % S00 = [SimGeom.H00plt(rngM,rngN), -SimGeom.G00plt(rngM,rngN)]*SimGeom.ZT + SimGeom.U00Tplt ...
    %     + SimGeom.G00plt(rngM,rngN)*SimGeom.QiT - SimGeom.H00plt(rngM,rngN)*SimGeom.UiT;
    % z00 = [SimGeom.H00plt(rngM,rngN), -SimGeom.G00plt(rngM,rngN)]*SimGeom.eT;       
    % SimSol.TL = S00*alphaT + z00;
    
    % uq = SimGeom.ZT*alphaT + SimGeom.eT;
    % N = length(uq)/2;
    % SimSol.TN = uq(1:N);
    % SimSol.Tknwn = [SimSol.TN; SimSol.TL];
    
    %--- SC Currents
    J_SC = struct('idx', {});
    for ii = 1:length(SimGeom.SC_tape)
        J_SC(ii).data = SimGeom.C3(SimGeom.SC_tape(ii).idx)/options.d_tape;
    end
    SimSol.J_SC = J_SC;

    %--- Domain Currents
    SimSol.jknwn = (SimGeom.Phi*alpha)./SimGeom.mu;
end