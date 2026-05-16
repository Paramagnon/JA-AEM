function [Jacobian] = returnJacobian(Geom,SimGeom,options,c3)
    % c3 vector stores the boundary conditions which are the J_i elements of
    % the superconductor. So c3 is zero if there is no superconductor element
    % between the boundaries. 
    % -------------------------------------------------------------------------
    Ec = options.Ec;
    n = options.n;
    Jc = options.Jc;    
    Jmin = Jc*1e-11;
    % Imin = 1e-12;

    d_tape = options.d_tape;       
    sigma_min = 1e5;
    EJ = @(j) min((Ec/Jc) * n * (max(abs(j),Jmin)/Jc).^(n-1),1/sigma_min);
     
    % -------------------------------------------------------------------------
    N = sum([Geom.Nbnd]);
    Nsc_tot = SimGeom.Nsc_tot; 
    % Nsc_tapes = SimGeom.Nsc_tapes;    
    Jacobian = SimGeom.JacobianInit;
    % C = zeros(Nsc_tapes,1);
    % -------------------------------------------------------------------------
    
    % -------------------------------------------------------------------------
    % Split jacobian in three contributions 
    % -------------------------------------------------------------------------       
    % J(Itot-Iset) 
    % Additional Lagrange multiplier per tape to ensure total current
    % follows the set current
    % Noffset = 0;
    % SC_it = 1;
    % for ii = 1:length(Geom)
    %     for nn = 1:length(Geom(ii).SC)
    %         SC_idx = Geom(ii).SC(nn).idx + Noffset;
    %         Atape = Geom(ii).SC(nn).dy*d_tape;
    %         Iset = currentSC(options,Geom(ii).SC(nn).pol);                       
    %         gamma = abs(Iset)/(abs(Iset) + Imin);
    % 
    %         J_curr = c3(SC_idx)/d_tape;            
    %         C(SC_it) = gamma*(Atape*J_curr-Iset);
    %         % disp(['Itot = ' num2str(Atape*sum(J_curr)) ', Iset = ' num2str(Iset)])
    %         % disp(['Atape = ' num2str(Atape) ', d_tape = ' num2str(d_tape)])
    % 
    %         SC_it = SC_it + 1;
    %     end
    %     Noffset = Noffset + Geom(ii).Nbnd;
    % end 
       
    % J(E_SC) - Only contributes on the diagonal elements where there
    % is a superconductor present, uses the E-J law of the SC 
    Noffset = 0;
    Nscoffset = 0;
    for ii = 1:length(Geom)
        for nn = 1:length(Geom(ii).SC)
            SC_idx = Geom(ii).SC(nn).idx;
            Nsc = length(SC_idx);
            for ll = 1:Nsc
                Ji = c3(Noffset + SC_idx(ll))/d_tape;
                Jacobian(Nscoffset+ll,Nscoffset+ll) = Jacobian(Nscoffset+ll,Nscoffset+ll) + EJ(Ji); 
            end   
            Nscoffset = Nscoffset + Nsc;
        end
        Noffset = Noffset + Geom(ii).Nbnd;
    end 
    
    %-------------------
    % Part below does not change, so no need to recalculate every time
    %------------------   
    % J(dA/dt)  
    % e = [H -G; C1 C2]^-1{0; C3}, so Jac(e) relates to rows in this matrix which for
    % simplicity is set as K = [H -G; C1 C2]^-1
    % invK = SimGeom.invK;
    % alpha(t+1) is related to C3 so Jac(alpha(t+1))    
    % AZ = SimGeom.AZ;
    %------------------
    Noffset = 0;
    Nscoffset = 0;
    for ii = 1:length(Geom)        
        for nn = 1:length(Geom(ii).SC)
            SCidxN = Geom(ii).SC(nn).idx + Noffset;
            SCidxJ = [1:length(SCidxN)] + Nscoffset;           
            %---------------------------------
            % A = Z \alpha + e
            % Jac(A) = Z Jac(\alpha) + Jac(e)
            %---------------------------------

            % --- Jac(\alpha) 
            % Using Ohm's law formulation:           
            Jacobian(SCidxJ,1:Nsc_tot) = Jacobian(SCidxJ,1:Nsc_tot) + options.d_tape*SimGeom.AZ(SCidxN,SimGeom.SCidx+N)/options.dt;
            % ---- Jac(e)
            % Using Ohm's law formulation:           
            % Jacobian(SCidxJ,1:Nsc_tot) = Jacobian(SCidxJ,1:Nsc_tot) + invK(SCidxN,SimGeom.SCidx+N)/options.dt;

            Nscoffset = Nscoffset + length(SCidxN);
        end
        Noffset = Noffset + Geom(ii).Nbnd;
    end
    Jacobian = mldivide(Jacobian,eye(size(Jacobian)));
end