function [boundary] = createRandomBoundaryAEM(p,bndr,options)
    boundary.x0 = p.x0;
    boundary.y0 = p.y0;
    
    boundary.domain_idx = p.domain_idx;
    boundary.boundary_idx = p.boundary_idx;
    boundary.element_idx = p.element_idx;
    
    boundary.xgeom = p.xgeom;
    boundary.ygeom = p.ygeom; 

    [boundary.Xm, boundary.Ym, boundary.bndr_idx, boundary.bndr_idxT, boundary.Nidx, boundary.c1, boundary.c2, boundary.c3, boundary.c1T, boundary.c2T, boundary.c3T] = returnGeomPointsAEM(boundary,bndr,options);
end