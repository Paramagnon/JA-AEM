function [JacobianInit, Nsc_tot, Nsc_tapes] = returnInitJacobian(Geom, SimGeom, options)
    d_tape = options.d_tape; 
    N = sum([Geom.Nbnd]);
    
    Nsc_tot = 0;
    Nsc_tapes = 0;
    for kk = 1:length(Geom)
        for nn = 1:length(Geom(kk).SC)
            Nsc_tapes = Nsc_tapes + 1;
            Nsc_tot = Nsc_tot + length(Geom(kk).SC(nn).idx);
        end
    end
    JacobianInit = zeros(Nsc_tot+Nsc_tapes,Nsc_tot+Nsc_tapes);
    
    %--- Some contributions to the Jacobian do not change so they are
    %initialized here

    % Additional Lagrange multiplier per tape to ensure total current
    % follows the set current
    Nscoffset = 0;
    Noffset = 0;
    SC_it = 1;
    for ii = 1:length(Geom)
        for nn = 1:length(Geom(ii).SC)
            SC_idx = Geom(ii).SC(nn).idx + Noffset;
            Nsc = length(SC_idx);
            Atape = Geom(ii).SC(nn).dy*d_tape;

            JacobianInit(Nsc_tot+SC_it, Nscoffset+1:Nscoffset+Nsc) = Atape';
            JacobianInit(Nscoffset+1:Nscoffset+Nsc,Nsc_tot+SC_it) = Atape;
           
            Nscoffset = Nscoffset + Nsc;
            SC_it = SC_it + 1;
        end
        Noffset = Noffset + Geom(ii).Nbnd;
    end 

    %--- Magnetic vector potential coupling
    invK = SimGeom.invK;        
    Noffset = 0;
    Nscoffset = 0;
    for ii = 1:length(Geom)        
        for nn = 1:length(Geom(ii).SC)
            SCidxN = Geom(ii).SC(nn).idx + Noffset;
            SCidxJ = [1:length(SCidxN)] + Nscoffset;           

            % Using Ohm's law formulation:           
            JacobianInit(SCidxJ,1:Nsc_tot) = JacobianInit(SCidxJ,1:Nsc_tot) + options.d_tape*invK(SCidxN,SimGeom.SCidx+N)/options.dt;

            Nscoffset = Nscoffset + length(SCidxN);
        end
        Noffset = Noffset + Geom(ii).Nbnd;
    end
    if nnz(JacobianInit)/numel(JacobianInit) < 0.1   
        JacobianInit = sparse(JacobianInit);
    end
end