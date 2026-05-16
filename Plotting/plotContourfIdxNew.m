function [idx, xknwn, yknwn, xp, yp] = plotContourfIdxNew(GeomN, GeomL, dlplot)
    xknwnN = zeros(1,sum([GeomN(1:end).Nbnd]));
    yknwnN = zeros(1,sum([GeomN(1:end).Nbnd]));
    xknwnL = zeros(1,sum([GeomL(1:end).Ldmn]));
    yknwnL = zeros(1,sum([GeomL(1:end).Ldmn]));

    Noffset = 0;
    Loffset = 0;
    for kk = 1:length(GeomL)
        Lpart = GeomL(kk).Ldmn;
        for ii = 1:length(GeomN(kk).boundary)
            Npart = length(GeomN(kk).boundary(ii).bndr_idx);

            xknwnN(Noffset+1:Noffset+Npart) = GeomN(kk).boundary(ii).Xn;
            yknwnN(Noffset+1:Noffset+Npart) = GeomN(kk).boundary(ii).Yn;
            Noffset = Noffset + Npart;
        end
        xknwnL(Loffset+1:Loffset+Lpart) = GeomL(kk).domain.Xn;
        yknwnL(Loffset+1:Loffset+Lpart) = GeomL(kk).domain.Yn;
        Loffset = Loffset + Lpart;
    end

    xmin = min(xknwnN);
    xmax = max(xknwnN);
    ymin = min(yknwnN);
    ymax = max(yknwnN);

    xi = linspace(xmin,xmax,round((xmax-xmin)*sqrt(dlplot),0));
    yi = linspace(ymin,ymax,round((ymax-ymin)*sqrt(dlplot),0));

    [Xi,Yi] = meshgrid(xi, yi);
    idx = zeros(size(Yi));

    xknwn = [xknwnN, xknwnL];
    yknwn = [yknwnN, yknwnL];
    xp = Xi;
    yp = Yi;
end