function [domain, Ldmn] = insertDomainPoints(Geom, options)
    % Input is one geom
    % dL = options.dL;
    dLx = options.dLx;
    dLy = options.dLy;

    rmin = options.rmin;
    % dL specifies a nodes density, with units [nodes/ area]
    
    % Find min and max x,y to create envelope encompassing the full
    % boundary
    xmin = Geom.boundary(1).Xm(1);
    xmax = Geom.boundary(1).Xm(1);
    ymin = Geom.boundary(1).Ym(1);
    ymax = Geom.boundary(1).Ym(1);

    for kk = 1:length(Geom.boundary)
        for ii = 1:length(Geom.boundary(kk).Xm)
            if Geom.boundary(kk).Xm(ii) > xmax
                xmax = Geom.boundary(kk).Xm(ii);
            elseif Geom.boundary(kk).Xm(ii) < xmin
                xmin = Geom.boundary(kk).Xm(ii);
            end
            if Geom.boundary(kk).Ym(ii) > ymax
                ymax = Geom.boundary(kk).Ym(ii);
            elseif Geom.boundary(kk).Ym(ii) < ymin
                ymin = Geom.boundary(kk).Ym(ii);
            end
        end
    end
    
    % Create square distribution of nodes based on dL in options, the mesh
    % density, equal square distribution

    % xi = linspace(xmin,xmax,round((xmax-xmin)*sqrt(dL),0));
    % yi = linspace(ymin,ymax,round((ymax-ymin)*sqrt(dL),0));

    xmin = xmin + 1.1*rmin;
    xmax = xmax - 1.1*rmin;
    ymin = ymin + 1.1*rmin;
    ymax = ymax - 1.1*rmin;
    xi = linspace(xmin,xmax,round((xmax-xmin)*dLx,0));
    yi = linspace(ymin,ymax,round((ymax-ymin)*dLy,0));
    
    [Xi,Yi] = meshgrid(xi,yi); 
    Xi = Xi(:)';
    Yi = Yi(:)';
    
    % Find nodes that are outside the boundaries, use the idea of ray
    % casting and seeing through how many edges the ray travels, if it is
    % uneven the point inside and vice verse. Draw line parallel to x-axis
    % from y point of the node, simplifies checks 

    idx = zeros(1,length(Xi));
    Boundary = Geom.boundary;

    for ii = 1:length(Yi)
        passes = 0;
        for kk = 1:length(Boundary)
            for nn = 1:length(Boundary(kk).Ym)-1
                Y1 = Boundary(kk).Ym(nn);
                Y2 = Boundary(kk).Ym(nn+1);                
                dY = Y2-Y1;
                if Yi(ii) > min([Y1,Y2]) && Yi(ii) <= max([Y1,Y2]) && dY ~= 0
                    X1 = Boundary(kk).Xm(nn);
                    X2 = Boundary(kk).Xm(nn+1);                  
                    Xint = (Yi(ii)-Y1)*(X2-X1)/(Y2-Y1) + X1;
                    if Xi(ii) > Xint
                        passes = passes+1;
                    end
                end           
            end
        end
        idx(ii) = mod(passes,2);
    end
    Xi = Xi(idx==1);
    Yi = Yi(idx==1);
    % Find nodes which do not satisfy the min distance condition to the
    % boundary
    
    idx = ones(1,length(Xi));
    
    for ii = 1:length(Yi)
        for kk = 1:length(Boundary)
            for nn = 1:length(Boundary(kk).Ym)-1
                Y1 = Boundary(kk).Ym(nn);
                Y2 = Boundary(kk).Ym(nn+1);
                X1 = Boundary(kk).Xm(nn);
                X2 = Boundary(kk).Xm(nn+1);

                A = Xi(ii) - X1;
                B = Yi(ii) - Y1;
                C = X2 - X1;
                D = Y2 - Y1;

                dist_case = (A*C + B*D)/(C*C + D*D);
                if dist_case < 0
                    xx = X1;
                    yy = Y1;
                elseif dist_case > 1
                    xx = X2;
                    yy = Y2;
                else
                    xx = X1 + dist_case*C;
                    yy = Y1 + dist_case*D;
                end

                dx = Xi(ii) - xx;
                dy = Yi(ii) - yy;
                r = sqrt(dx*dx + dy*dy);
                if r < rmin
                    idx(ii) = 0;                
                end           
            end
        end
    end   
    domain.Xn = Xi(idx==1);
    domain.Yn = Yi(idx==1);
    domain.rmin = rmin;
    Ldmn = length(domain.Xn);
end
    