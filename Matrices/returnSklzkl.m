function [S00,Sx0,Sy0,z00,zx0,zy0] = returnSklzkl(SimGeom)
    % [S00,Sx0,Sy0,Sxx,Syy,z00,zx0,zy0,zxx,zyy]

    % Skl = [Hkl -Gkl] Z + (Ui + Gkl Q - Hkl U)
    S00 = [SimGeom.H00, -SimGeom.G00]*SimGeom.Z + SimGeom.U00 + SimGeom.G00*SimGeom.Qi - SimGeom.H00*SimGeom.Ui;
    Sx0 = [SimGeom.Hx0, -SimGeom.Gx0]*SimGeom.Z + SimGeom.Ux0 + SimGeom.Gx0*SimGeom.Qi - SimGeom.Hx0*SimGeom.Ui;
    Sy0 = [SimGeom.Hy0, -SimGeom.Gy0]*SimGeom.Z + SimGeom.Uy0 + SimGeom.Gy0*SimGeom.Qi - SimGeom.Hy0*SimGeom.Ui;
    % Sxx = [Geom.Hxx, -Geom.Gxx]*Geom.Z + Geom.Uxx + Geom.Gxx*Geom.Qi - Geom.Hxx*Geom.Ui;
    % Syy = [Geom.Hyy, -Geom.Gyy]*Geom.Z + Geom.Uyy + Geom.Gyy*Geom.Qi - Geom.Hyy*Geom.Ui;

    % zkl = [Hkl -Gkl] e
    z00 = [SimGeom.H00, -SimGeom.G00]*SimGeom.e;
    zx0 = [SimGeom.Hx0, -SimGeom.Gx0]*SimGeom.e;
    zy0 = [SimGeom.Hy0, -SimGeom.Gy0]*SimGeom.e;
    % zxx = [Geom.Hxx, -Geom.Gxx]*Geom.e;
    % zyy = [Geom.Hyy, -Geom.Gyy]*Geom.e;

end