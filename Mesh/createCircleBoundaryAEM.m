function [boundary] = createCircleBoundaryAEM(p,bndr,options,cw)
    boundary.x0 = p.x0;
    boundary.y0 = p.y0;
    
    boundary.domain_idx = p.domain_idx;
    boundary.boundary_idx = p.boundary_idx;
    boundary.element_idx = p.element_idx;       

    [boundary.Xm, boundary.Ym, boundary.bndr_idx, boundary.bndr_idxT, boundary.Nidx, boundary.c1, boundary.c2, boundary.c3, boundary.c1T, boundary.c2T, boundary.c3T] =  returnMeshCircleAEM(p,bndr,options,cw);
    boundary.xgeom = [boundary.Xm(1), boundary.Xm(2)];
    boundary.ygeom = [boundary.Ym(1), boundary.Ym(2)];
end