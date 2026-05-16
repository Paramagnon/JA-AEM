function [SimGeom] = runSteadyState(Geom, SimGeom, options)
       
    [SimGeom.SCidx, SimGeom.SC_tape] = returnVectorSCidx(Geom);
    
    [c1,c2,SimGeom.C3] = returnVectorc1c2c3(Geom);
    SimGeom.C1 = SimGeom.C1idx.*c1';
    SimGeom.C2 = SimGeom.C2idx.*c2';

    [SimGeom.Z,SimGeom.e] = returnMatrixZe(SimGeom);
    [SimGeom.S00,SimGeom.Sx0,SimGeom.Sy0,SimGeom.z00,SimGeom.zx0,SimGeom.zy0] = returnSklzkl(SimGeom);
        
    SimGeom.f = zeros(length(SimGeom.S00),1);
    alpha = (SimGeom.Uxx + SimGeom.Uyy)\SimGeom.f;

    options.t = 0;
    SimGeom.SimSol = returnSimSol(SimGeom,alpha,options);   
    
end