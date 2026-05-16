function [Geom, plt] = plotContourfU(Geom, SimSol, options, reinit)    
    dlplot = options.dlplot; 
    crcl = options.crcl;
    
    if reinit
        [plt.idx, plt.xknwn, plt.yknwn, plt.xp, plt.yp] = plotContourfIdxNew(Geom, dlplot);
    else 
        plt = options.plt;
    end
    
    warning('off','MATLAB:griddata:DuplicateDataPoints')
    plotgridu = griddata(plt.xknwn, plt.yknwn, SimSol.uknwn, plt.xp, plt.yp,"cubic");    
    contourf(plt.xp, plt.yp, plotgridu, 20);
    
    if crcl
        hold on
        if length(Geom(end).boundary) > 1
            for ii = 2:length(Geom(end).boundary)
                r = Geom(end).boundary(ii).xgeom(1)-Geom(end).boundary(ii).x0;
                circles(Geom(end).boundary(ii).x0,Geom(end).boundary(ii).y0,r,"facecolor","black")
            end
        end
        hold off
    end
end