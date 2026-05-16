function [e,z00,zx0,zy0] = updateC3(SimGeom)

    % e = [Geom.H, - Geom.G; Geom.C1, Geom.C2] \ [Geom.p; Geom.C3];
    e = SimGeom.invK*[SimGeom.p; SimGeom.C3];

    z00 = [SimGeom.H00, -SimGeom.G00]*e;
    zx0 = [SimGeom.Hx0, -SimGeom.Gx0]*e;
    zy0 = [SimGeom.Hy0, -SimGeom.Gy0]*e;

end