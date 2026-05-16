function [plt] = plotContourfThermal(GeomN, GeomL, SimSol, options, reinit)      
    dlplot = options.dlplot; 
    % crcl = options.crcl;
    % N = sum([GeomN.Nbnd]);

    if reinit 
        [plt.idx, plt.xknwn, plt.yknwn, plt.xp, plt.yp] = plotContourfIdxNew(GeomN, GeomL, dlplot);
    else 
        plt = options.plt;
    end
            
    warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
    % if options.pltL
    %     plotgridu = griddata(plt.xknwn(N+1:end), plt.yknwn(N+1:end), SimSol.jknwn, plt.xp, plt.yp,"cubic");
    %     image([min(plt.xp),max(plt.xp)], [min(plt.yp),max(plt.yp)], plotgridu,'CDataMapping','scaled')
    %     % image(plt.xp, plt.yp, plotgridu,"EdgeColor","none","FaceColor","interp")
    %     hold on
    %     contour(plt.xp, plt.yp, plotgridu, 20,"EdgeColor","k");
    %     hold off
    %     ax = gca;
    %     ax.YDir = 'normal';
    % else
    F1 = scatteredInterpolant(plt.xknwn',plt.yknwn',SimSol.Tknwn,'linear','linear');
    plotgridu = F1(plt.xp,plt.yp);
    % image([min(plt.xp),max(plt.xp)], [min(plt.yp),max(plt.yp)], plotgridu,'CDataMapping','scaled')
    hold on
    contourf(plt.xp, plt.yp, plotgridu, 20,"EdgeColor","k");
    hold off
    ax = gca;
    ax.YDir = 'normal';

    % hold on   
    % if options.pltlines
    %     for nn = 1:3
    %         for ii = 1:length(GeomN(nn).boundary(1).xgeom)-1
    %             x0 = GeomN(nn).boundary(1).xgeom(ii);
    %             x1 = GeomN(nn).boundary(1).xgeom(ii+1);
    %             y0 = GeomN(nn).boundary(1).ygeom(ii);
    %             y1 = GeomN(nn).boundary(1).ygeom(ii+1);
    %             plot([x0,x1],[y0,y1],"k--", "LineWidth",2)
    %         end
    %     end
    % end
    % hold off
    % if crcl
    %     hold on
    %     if length(Geom(end).boundary) > 1
    %         for ii = 2:length(Geom(end).boundary)
    %             r = Geom(end).boundary(ii).xgeom(1)-Geom(end).boundary(ii).x0;
    %             circles(Geom(end).boundary(ii).x0,Geom(end).boundary(ii).y0,r,"facecolor","black")
    %         end
    %     end
    %     hold off
    % end
end