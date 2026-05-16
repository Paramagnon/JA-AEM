function [xgeom, ygeom] = returnSquarePoints(p, cw)   
    % Define corners as x0, x1, y0, y1
    % Returns points for inner domain
    % cw = 1 is clockwise domain, cw = 0 is anti clockwise

    xgeom = [p.x0,p.x1,p.x1,p.x0,p.x0];
    ygeom = [p.y0,p.y0,p.y1,p.y1,p.y0]; 

    if cw 
        xgeom = xgeom(end:-1:1);
        ygeom = ygeom(end:-1:1);
    end
end