function [C] = returnC(Geom,SimGeom,options,c3)
    % c3 vector stores the boundary conditions which are the J_i elements of
    % the superconductor. So c3 is zero if there is no superconductor element
    % between the boundaries. 
    % -------------------------------------------------------------------------
    Imin = 1e-12;
    d_tape = options.d_tape;       
    % -------------------------------------------------------------------------
    Nsc_tapes = SimGeom.Nsc_tapes;    
    C = zeros(Nsc_tapes,1);
    % ------------------------------------------------------------------------- 
    % J(Itot-Iset) 
    % Additional Lagrange multiplier per tape to ensure total current
    % follows the set current
    Noffset = 0;
    SC_it = 1;
    for ii = 1:length(Geom)
        for nn = 1:length(Geom(ii).SC)
            SC_idx = Geom(ii).SC(nn).idx + Noffset;
            Atape = Geom(ii).SC(nn).dy*d_tape;
            Iset = currentSC(options,Geom(ii).SC(nn).pol);                       
            gamma = abs(Iset)/(abs(Iset) + Imin);

            J_curr = c3(SC_idx)/d_tape;            
            C(SC_it) = gamma*(Atape*J_curr-Iset);
            % disp(['Itot = ' num2str(Atape*sum(J_curr)) ', Iset = ' num2str(Iset)])
            % disp(['Atape = ' num2str(Atape) ', d_tape = ' num2str(d_tape)])
            SC_it = SC_it + 1;
        end
        Noffset = Noffset + Geom(ii).Nbnd;
    end 
end