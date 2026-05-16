function [S00,z00] = returnSklzklThermal(SimGeom)
    % [S00,Sx0,Sy0,Sxx,Syy,z00,zx0,zy0,zxx,zyy]

    % Skl = [Hkl -Gkl] Z + (Ui + Gkl Q - Hkl U)
    S00 = [SimGeom.H00T, -SimGeom.G00T]*SimGeom.ZT + SimGeom.U00T + SimGeom.G00T*SimGeom.QiT - SimGeom.H00T*SimGeom.UiT;
    % Sx0 = [Geom.Hx0, -Geom.Gx0]*Geom.Z + Geom.Ux0 + Geom.Gx0*Geom.Qi - Geom.Hx0*Geom.Ui;
    % Sy0 = [Geom.Hy0, -Geom.Gy0]*Geom.Z + Geom.Uy0 + Geom.Gy0*Geom.Qi - Geom.Hy0*Geom.Ui;
    % Sxx = [Geom.Hxx, -Geom.Gxx]*Geom.Z + Geom.Uxx + Geom.Gxx*Geom.Qi - Geom.Hxx*Geom.Ui;
    % Syy = [Geom.Hyy, -Geom.Gyy]*Geom.Z + Geom.Uyy + Geom.Gyy*Geom.Qi - Geom.Hyy*Geom.Ui;

    % zkl = [Hkl -Gkl] e
    z00 = [SimGeom.H00T, -SimGeom.G00T]*SimGeom.eT;
    % zx0 = [Geom.Hx0, -Geom.Gx0]*Geom.e;
    % zy0 = [Geom.Hy0, -Geom.Gy0]*Geom.e;
    % zxx = [Geom.Hxx, -Geom.Gxx]*Geom.e;
    % zyy = [Geom.Hyy, -Geom.Gyy]*Geom.e;

end