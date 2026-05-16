function [meshPointsX, meshPointsY] = returnMeshLinearAEM(x0,y0,x1,y1,N)     
    meshPointsX = linspace(x0,x1,N);
    meshPointsY = linspace(y0,y1,N);
end

