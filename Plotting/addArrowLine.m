function addArrowLine(xgeom,ygeom)
    xp = (xgeom(1) + xgeom(2))/2;
    yp = (ygeom(1) + ygeom(2))/2;
    dx = xgeom(2) - xgeom(1);
    dy = ygeom(2) - ygeom(1);

    quiver(xp,yp,dx,dy,0.25,"Color","k","LineWidth",4,"AutoScale","on","MaxHeadSize",12)
end