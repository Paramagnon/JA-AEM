function [alpha,SimGeom,PQ] = runMagnetic(SimGeom, Geom, alpha_prev, options)
    % Find SC current distribution which matches the next time step and
    % satisfies a pre-defined accuracy value   
    PQ = 0;
    % -------------------------------------------------------------------------
    Ec = options.Ec;
    n = options.n;
    Jc = options.Jc;
    d_tape = options.d_tape;
    omega = 1; % Damping for updating current densities
    kjacup = 3; % Update Jacobian only every so many iterations
    % -------------------------------------------------------------------------     
    %---------------- Matrices ---------------------
    dtU = options.dt*SimGeom.UxxUyy;
    dtf = options.dt*(SimGeom.f); % Misses f_prev term, not needed if Fourier Series is used
    C0 = SimGeom.mu.*SimGeom.sigmaA;
    CS = C0.*SimGeom.S00;
    if options.Tupdate
        if options.Euler
            % Backward Euler:
            invdtUCS = mldivide((dtU - CS)*1e6,eye(size(dtU)))*1e6;
        else
            % Crank–Nicolson:
            invdtUCS = mldivide((dtU - 2*CS)*1e6,eye(size(dtU)))*1e6;
        end
        SimGeom.AZ = SimGeom.Z*invdtUCS*(C0.*[SimGeom.H00, -SimGeom.G00]*SimGeom.invK);   
    else
        SimGeom.AZ = SimGeom.Z*SimGeom.invdtUCS*(C0.*[SimGeom.H00, -SimGeom.G00]*SimGeom.invK);   
    end
    % -------------------------------------------------------------------------
    sigma_min = 1e5;
    EJ = @(j) sign(j).*min((Ec*(abs(j)/Jc).^n),Jc/sigma_min); 
    %---------------- Number of SC tapes -----------
    N = sum([Geom.Nbnd]);    
    Nsc_tot = SimGeom.Nsc_tot; 
    if Nsc_tot == 0
        if options.Tupdate
            if options.Euler
                % Backward Euler:
                alpha = invdtUCS*(dtf - CS*alpha_prev);
            else
                % Crank–Nicolson:
                alpha = invdtUCS*(dtf - (2*CS + dtU)*alpha_prev);
            end
        else
            if options.Euler
                % Backward Euler:
                alpha = SimGeom.invdtUCS*(dtf - CS*alpha_prev);
            else               
                % Crank–Nicolson:
                alpha = SimGeom.invdtUCS*(dtf - (2*CS + dtU)*alpha_prev);
            end
        end
        return
    end 
    %---------------- Initialize variables ---------
    p_prev = SimGeom.p;    
    [SimGeom.puN, SimGeom.puL, SimGeom.pqN] = boundaryFourierCurrent(SimGeom,Geom,options,options.t);
    SimGeom.p = SimGeom.H*SimGeom.puN - SimGeom.G*SimGeom.pqN;

    g = zeros(Nsc_tot,1);
    C3_prev = SimGeom.C3;
    C3_curr = C3_prev;
    % PQ = zeros(size(C3_curr));                              % Heat source on the boundary
    C3_full = [SimGeom.p; SimGeom.C3];
    dC3 = zeros(2*N,1);
    dC3(1:N) = SimGeom.p - p_prev;

    J_curr = C3_curr(SimGeom.SCidx)/d_tape;   
    %---------------- Previous values --------------
    uq = SimGeom.Z*alpha_prev + SimGeom.e;
    A_prev = uq(1:N);
    %---------------- First J guess ----------------
    [invJacobian] = returnJacobian(Geom,SimGeom,options,C3_curr);
    [C] = returnC(Geom,SimGeom,options,C3_curr);
    dJ = invJacobian*[g;C];
    if (sum(isnan(dJ))~= 0)
        w = SimGeom.JacobianInit(Nsc_tot+1:end,1:Nsc_tot);
        J_curr = (sum(w.*(C./sum(w,2)),1))./sum(w,1);
        % J_curr = ((C.^-1)'*SimGeom.JacobianInit(Nsc_tot+1:end,1:Nsc_tot)).^-1; %
    else
        J_curr = J_curr - omega*dJ(1:Nsc_tot);
    end
    C3_curr(SimGeom.SCidx) = J_curr*d_tape;
    
    for kk = 2:options.maxit       
        dC3(N+1:2*N) = C3_curr-C3_prev;        
        dz = ([SimGeom.H00, -SimGeom.G00]*SimGeom.invK)*dC3;
        if options.Tupdate
            if options.Euler
                % Backward Euler:
                alpha = invdtUCS*(dtf + C0.*dz - CS*alpha_prev);
            else
                % Crank–Nicolson: 
                alpha = invdtUCS*(dtf + 2*C0.*dz - (2*CS + dtU)*alpha_prev);
            end
        else
            if options.Euler
                % Backward Euler:
                alpha = SimGeom.invdtUCS*(dtf + C0.*dz - CS*alpha_prev);
            else
                % Crank–Nicolson: 
                alpha = SimGeom.invdtUCS*(dtf + 2*C0.*dz - (2*CS + dtU)*alpha_prev);           
            end
        end

        C3_full(N+1:2*N) = C3_curr;
        e = SimGeom.invK*C3_full;       
        uq = SimGeom.Z*alpha + e;
        Nuq = length(uq)/2;
        A = uq(1:Nuq);

        Noffset = 0;
        Nscoffset = 0;
        for ii = 1:length(Geom)
            for nn = 1:length(Geom(ii).SC)
                SCidxN = Geom(ii).SC(nn).idx + Noffset;
                Nsc = length(SCidxN);

                Ayy = (A(SCidxN) - A_prev(SCidxN))/options.dt;
                Jsc = C3_curr(SCidxN)/d_tape;                               
                
                gSC = EJ(Jsc) + Ayy;
                
                g(Nscoffset+1:Nscoffset+Nsc) = gSC - mean(gSC);

                Nscoffset = Nscoffset + Nsc;
            end
            Noffset = Noffset + Geom(ii).Nbnd;
        end
        if sum(abs(g) < options.maxerror) == length(g)
            C3_curr(SimGeom.SCidx) = J_curr*d_tape;           
            SimGeom.C3 = C3_curr;

            % PQ(SimGeom.SCidx) = EJ(J_curr).*J_curr;                

            [SimGeom.e,SimGeom.z00,SimGeom.zx0,SimGeom.zy0] = updateC3(SimGeom);
            return
        else                       
            if mod(kk,kjacup) == 0
                [invJacobian] = returnJacobian(Geom,SimGeom,options,C3_curr);
            end
            [C] = returnC(Geom,SimGeom,options,C3_curr);
            dJ = invJacobian*[g;C];
            if (sum(isnan(dJ))~= 0)
                w = SimGeom.JacobianInit(Nsc_tot+1:end,1:Nsc_tot);
                J_curr = (sum(w.*(C./sum(w,2)),1))./sum(w,1);
                % J_curr = ((C.^-1)'*SimGeom.JacobianInit(Nsc_tot+1:end,1:Nsc_tot)).^-1; %
            else
                J_curr = J_curr - omega*dJ(1:Nsc_tot);
            end
            C3_curr(SimGeom.SCidx) = J_curr*d_tape;
            % if options.debug
            %     disp(['Newton-Raphson iteration: ' num2str(kk)])
            %     disp(g')
            % end 
        end
    end
    error("Max Netwon-Raphson itterations reached without convergence");
end
