function [Geom] = returnNodesAEM(Geom)   
    for kk = 1:length(Geom)
        for ii = 1:length(Geom(kk).boundary)
            X = Geom(kk).boundary(ii).Xm;
            Y = Geom(kk).boundary(ii).Ym;
    
            Geom(kk).boundary(ii).Xn = (X(1:end-1)+X(2:end))/2; % Xn is one shorter than X, X double counts x0, node coordinates
            Geom(kk).boundary(ii).Yn = (Y(1:end-1)+Y(2:end))/2; % Yn is one shorter than Y, Y double counts y0, node coordinates
        end
    end
end