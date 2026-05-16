function [Z,e] = returnMatrixZeThermal(SimGeom)
    Trng = SimGeom.Rngs.Nrng;  
    Z = [SimGeom.H(Trng,Trng), - SimGeom.G(Trng,Trng); SimGeom.C1T, SimGeom.C2T] \ [SimGeom.TT; zeros(size(SimGeom.TT))];
    e = [SimGeom.H(Trng,Trng), - SimGeom.G(Trng,Trng); SimGeom.C1T, SimGeom.C2T] \ [zeros(size(SimGeom.C3T)); SimGeom.C3T];

end