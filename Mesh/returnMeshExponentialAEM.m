function [meshPointsX, meshPointsY] = returnMeshExponentialAEM(x0,y0,x1,y1,dmin,gamma)
    flipx = 0;
    flipy = 0;
    if x0 > x1
        temp = x1;
        x1 = x0;
        x0 = temp;        
        flipx = 1;
    end
    if y0 > y1
        temp = y1;
        y1 = y0;
        y0 = temp;
        flipy = 1;
    end
    
    if x0 == x1     
        y2 = (y1+y0)/2;              
        
        dy = (dmin/(y2-y0))*(gamma*3-1);
        dy = log(1+dy) - log(1);
        
        meshPointsYint = (1/(gamma*3-1))*(exp(log(1):dy:log(gamma*3))-1)*(y2-y0)+y0;
        if (y2 - meshPointsYint(end) < meshPointsYint(end) - meshPointsYint(end-1))
            meshPointsY = [meshPointsYint(1:end-1), y2, 2*y2-meshPointsYint(end-1:-1:1)]; 
        else
            meshPointsY = [meshPointsYint, y2, 2*y2-meshPointsYint(end:-1:1)];
        end       
        meshPointsX = ones(size(meshPointsY))*x0;
    else       
        x2 = (x1+x0)/2;    
        alpha = (y1-y0)/(x1-x0);
    
        f_a = (1+alpha^2);
        f_b = -2*x0*(1+alpha^2);
        f_c = (1+alpha^2)*x0^2-dmin^2;        
        xb = (-f_b + sqrt(f_b^2 - 4*f_a*f_c))/(2*f_a);
                      
        dx = ((xb-x0)/(x2-x0))*(gamma*3-1);
        dx = log(1+dx) - log(1);
        
        meshPointsXint = (1/(gamma*3-1))*(exp(log(1):dx:log(gamma*3))-1)*(x2-x0)+x0;        
        if (x2 - meshPointsXint(end) < meshPointsXint(end) - meshPointsXint(end-1))
            meshPointsX = [meshPointsXint(1:end-1), x2, 2*x2-meshPointsXint(end-1:-1:1)]; 
        else
            meshPointsX = [meshPointsXint, x2, 2*x2-meshPointsXint(end:-1:1)]; 
        end               
        meshPointsY = ((y1-y0)/(x1-x0))*meshPointsX + y0 - ((y1-y0)/(x1-x0))*x0;       
    end    
    if flipx
        meshPointsX = meshPointsX(end:-1:1);
    end
    if flipy
        meshPointsY = meshPointsY(end:-1:1);
    end
end

