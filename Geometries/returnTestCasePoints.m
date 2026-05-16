function [xgeom, ygeom] = returnTestCasePoints(p)
    % Define hb, wb, hd, wd
    % Returns points for inner domain
    x0 = p.x0;
    y0 = p.y0;
    hb = p.hb;
    wb = p.wb;
    hd = p.hd;
    wd = p.wd;
       
    xgeom = [x0,x0,x0+wb,x0+wb,x0+wb+wd,x0+wb+wd,x0+2*wb+wd,x0+2*wb+wd,x0+2*wb+2*wd, ...
        x0+2*wb+2*wd,x0+3*wb+2*wd,x0+3*wb+2*wd,x0];
    ygeom = [y0,y0+hb,y0+hb,y0+hb-hd,y0+hb-hd,y0+hb,y0+hb,y0+hb-hd,y0+hb-hd,y0+hb, ...
        y0+hb,y0,y0]; 
    xgeom = xgeom(end:-1:1);
    ygeom = ygeom(end:-1:1);
end