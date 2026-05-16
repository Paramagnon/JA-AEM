function [Z,e] = returnMatrixZe(Geom)

    Z = [Geom.H, - Geom.G; Geom.C1, Geom.C2] \ [Geom.T; zeros(size(Geom.T))];
    e = [Geom.H, - Geom.G; Geom.C1, Geom.C2] \ [zeros(size(Geom.C3)); Geom.C3];

end