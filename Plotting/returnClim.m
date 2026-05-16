function [Cmin, Cmax] = returnClim(SimSol,options)
    if options.pltL
        Cmin = min(SimSol(1).jknwn);
        Cmax = max(SimSol(1).jknwn);
    
        if length(SimSol) > 1
            for tt = 2:length(SimSol)
                Cmin = min([SimSol(tt).jknwn; Cmin]);
                Cmax = max([SimSol(tt).jknwn; Cmax]);
            end
        end
        dc = Cmax - Cmin;
        if dc == 0
            % dc = 0.1;
            Cmin = -1;
            Cmax = 1;
        end
    else        
        Cmin = min([SimSol(1).uN; SimSol(1).uL]);
        Cmax = max([SimSol(1).uN; SimSol(1).uL]);
    
        if length(SimSol) > 1
            for tt = 2:length(SimSol)
                Cmin = min([SimSol(tt).uN; SimSol(tt).uL; Cmin]);
                Cmax = max([SimSol(tt).uN; SimSol(tt).uL; Cmax]);
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
end