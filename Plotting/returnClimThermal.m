function [Cmin, Cmax] = returnClimThermal(SimSol,options)
       
    Cmin = min([SimSol(1).TN; SimSol(1).TL]);
    Cmax = max([SimSol(1).TN; SimSol(1).TL]);
    
    if length(SimSol) > 1
        for tt = 2:length(SimSol)
            Cmin = min([SimSol(tt).TN; SimSol(tt).TL; Cmin]);
            Cmax = max([SimSol(tt).TN; SimSol(tt).TL; Cmax]);
        end
    end
    dc = Cmax - Cmin;
    if dc == 0
        % dc = 0.1;
        Cmin = -1;
        Cmax = 1;
    end
    % Cmin = Cmin - 0.2*dc;
    % Cmax = Cmax + 0.2*Cmax;

end