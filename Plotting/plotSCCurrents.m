function plotSCCurrents(SimGeom,options,solidx)       
    fileSC = 'SCcurrents.gif';
    
    SCtot = length(SimGeom.SimSol(1).J_SC);
    idx1 = 1;
    idx2 = 5;
    idx3 = 6;
    idx4 = 10;
    % idx2 = SCtot/8;
    % idx3 = SCtot/8 + 1;
    % idx4 = SCtot/4;
    
    figure("Position",[70,50,500,430],"Units","pixels")
    % figure
    % axis square
    hold on
    xlim([0,length(SimGeom.SimSol(1).J_SC(idx1).data)*options.dy])
    SC1 = stairs(nan,nan,"LineWidth",2,"Color","r","LineStyle","-");
    if SCtot > 8 
        SC2 = stairs(nan,nan,"LineWidth",2,"Color","y","LineStyle","-");
        SC3 = stairs(nan,nan,"LineWidth",2,"Color","g","LineStyle","-");
    end
    SC4 = stairs(nan,nan,"LineWidth",2,"Color","b","LineStyle","-");
    hold off

    if solidx == 0
        SC1.XData = (0.5:length(SimGeom.SimSol(1).J_SC(idx1).data)-0.5)*options.dy;
        SC1.YData = SimGeom.SimSol(1).J_SC(idx1).data/options.Jc;       
        if SCtot > 8
            SC2.XData = (0.5:length(SimGeom.SimSol(1).J_SC(idx2).data)-0.5)*options.dy;
            SC2.YData = SimGeom.SimSol(1).J_SC(idx2).data/options.Jc;
            SC3.XData = (0.5:length(SimGeom.SimSol(1).J_SC(idx3).data)-0.5)*options.dy;
            SC3.YData = SimGeom.SimSol(1).J_SC(idx3).data/options.Jc;
        end
        SC4.XData = (0.5:length(SimGeom.SimSol(1).J_SC(idx4).data)-0.5)*options.dy;
        SC4.YData = SimGeom.SimSol(1).J_SC(idx4).data/options.Jc;
            
        grid on
        xlabel("Tape width [m]","FontSize",12,"Interpreter","latex")
        ylabel("Current density J/Jc [-]","FontSize",12,"Interpreter","latex")
        if SCtot > 8 
            legend([SC1,SC2,SC3,SC4],['Tape 1: Iset = ' num2str(-options.Iset) ' A'],['Tape ' num2str(idx2) ': Iset = ' num2str(-options.Iset) ' A'],... 
                ['Tape ' num2str(idx3) ': Iset = ' num2str(options.Iset) ' A'],['Tape ' num2str(idx4) ': Iset = ' num2str(options.Iset) ' A'],"location","north","numcolumns",2)    
        else
            legend([SC1,SC4],['Tape 1: Iset = ' num2str(-options.Iset) ' A'],['Tape 2: Iset = ' num2str(options.Iset) ' A'],"location","north")
        end
        ylim([-1.3 1.7])
        title(['t = ' num2str(0*options.dt) ' s'])
        exportgraphics(gcf,fileSC);
        
        for tt = 2:length(SimGeom.SimSol)                    
            SC1.YData = SimGeom.SimSol(tt).J_SC(idx1).data/options.Jc;       
            if SCtot > 8
                SC2.YData = SimGeom.SimSol(tt).J_SC(idx2).data/options.Jc;
                SC3.YData = SimGeom.SimSol(tt).J_SC(idx3).data/options.Jc;       
            end
            SC4.YData = SimGeom.SimSol(tt).J_SC(idx4).data/options.Jc;
            exportgraphics(gcf,fileSC,'Append',true);        
            title(['t = ' num2str((tt-1)*options.dt) ' s'])
        end   
    
    else
        SC1.XData = (0.5:length(SimGeom.SimSol(solidx).J_SC(idx1).data)-0.5)*options.dy;
        SC1.YData = SimGeom.SimSol(solidx).J_SC(idx1).data/options.Jc;       
        if SCtot > 8
            SC2.XData = (0.5:length(SimGeom.SimSol(solidx).J_SC(idx2).data)-0.5)*options.dy;
            SC2.YData = SimGeom.SimSol(solidx).J_SC(idx2).data/options.Jc;
            SC3.XData = (0.5:length(SimGeom.SimSol(solidx).J_SC(idx3).data)-0.5)*options.dy;
            SC3.YData = SimGeom.SimSol(solidx).J_SC(idx3).data/options.Jc;
        end
        SC4.XData = (0.5:length(SimGeom.SimSol(solidx).J_SC(idx4).data)-0.5)*options.dy;
        SC4.YData = SimGeom.SimSol(solidx).J_SC(idx4).data/options.Jc;
        
        grid on
        xlabel("Tape width [m]","FontSize",12,"Interpreter","latex")
        ylabel("Current density J/Jc [-]","FontSize",12,"Interpreter","latex")
        if SCtot > 8 
            legend([SC1,SC2,SC3,SC4],['Tape 1: Iset = ' num2str(-options.Iset) ' A'],['Tape ' num2str(idx2) ': Iset = ' num2str(-options.Iset) ' A'],... 
                ['Tape ' num2str(idx3) ': Iset = ' num2str(options.Iset) ' A'],['Tape ' num2str(idx4) ': Iset = ' num2str(options.Iset) ' A'],"location","north","numcolumns",2)    
        else
            legend([SC1,SC4],['Tape 1: Iset = ' num2str(-options.Iset) ' A'],['Tape 2: Iset = ' num2str(options.Iset) ' A'],"location","north")
        end
        ylim([-1.3 1.7])
        title(['t = ' num2str((solidx-1)*options.dt) ' s'])

    end
end
       
    

    
