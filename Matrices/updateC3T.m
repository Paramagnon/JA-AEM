function [e,z00] = updateC3T(SimGeom)

    % e = [Geom.H, - Geom.G; Geom.C1, Geom.C2] \ [Geom.p; Geom.C3];
    e = SimGeom.invKT*[zeros(size(SimGeom.C3T)); SimGeom.C3T];
    z00 = [SimGeom.H00T, -SimGeom.G00T]*e;   

end