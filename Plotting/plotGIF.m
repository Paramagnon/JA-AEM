function plotGIF(Geom, SimGeom, options)
    t = 0:options.dt:options.tmax;
    [Cmin, Cmax] = returnClim(SimGeom.SimSol,options);

    % figure
    figure("Position",[1,1,828,499],"Units","pixels")

    if options.debug
        disp(['Plot timestep: ' num2str(1)])
    end
    [options.plt] = plotContourf(Geom, SimGeom.PlotGeom, SimGeom.SimSol(1), options, 1);
    colorMapAltair()
    c = colorbar;
    if options.pltL
        c.Label.String = 'Current density [A/m$^2$]';
    else
        c.Label.String = 'Magnitude [Wb/m]';
    end
    c.Label.FontSize = 12;
    c.Label.Interpreter = "latex";
    clim([Cmin Cmax])

    % axis square
    grid on
    xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
    ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
    exportgraphics(gcf,options.file);

    for tt = 2:length(t)
        if options.debug
            disp(['Plot timestep: ' num2str(tt)])
        end

        [~] = plotContourf(Geom, SimGeom.PlotGeom, SimGeom.SimSol(tt), options, 0);

        colorMapAltair()
        c = colorbar;
        if options.pltL
            c.Label.String = 'Current density [A/m$^2$]';
        else
            c.Label.String = 'Magnitude [Wb/m]';
        end
        c.Label.FontSize = 12;
        c.Label.Interpreter = "latex";
        clim([Cmin Cmax])

        % axis square
        grid on
        xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
        ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
        exportgraphics(gcf,options.file,'Append',true);
    end
end