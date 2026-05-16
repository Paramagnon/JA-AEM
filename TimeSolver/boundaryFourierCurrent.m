function [uN,uL,qN] = boundaryFourierCurrent(SimGeom,Geom,options,t)
    par = options.par;
    Jref = -4*pi*1e-7*currentAC(options, par);
    % Jref = currentAC(options, par);   
    k_n = @(n) n*pi/par.tau;

    uN = zeros(size(SimGeom.C3));
    uL = zeros(size(SimGeom.G00plt,1),1);
    qN = zeros(size(SimGeom.C3));

    %--- Motion profile
    theta = 2*pi*par.F*t;
    tau_m = 0.03;
    x = 4*tau_m*sin(theta)+4*tau_m;
    theta = pi*x/tau_m;
    %--- Coil currents
    Ja = cos(theta)*Jref;
    Jb = cos(theta-2*pi/3)*Jref;
    Jc = cos(theta+2*pi/3)*Jref;

    %--- Actual values
    Nn = [1:par.NFour]';  
    k_nn = k_n(Nn);
    k_nn_1 = 1./k_nn';
    k_nn_2 = k_nn_1.*k_nn_1;

    theta = 0; %---> Turn off movement, if you want pulsing coils
    costheta = cos(k_nn'*theta);
    sintheta = sin(k_nn'*theta);
    
    an = SimGeom.an;
    bn = SimGeom.bn;
    
    an_nn = Ja*(an(1,:).*costheta + bn(1,:).*sintheta) + Jb*(an(2,:).*costheta + bn(2,:).*sintheta) + Jc*(an(3,:).*costheta + bn(3,:).*sintheta);
    bn_nn = Ja*(bn(1,:).*costheta - an(1,:).*sintheta) + Jb*(bn(2,:).*costheta - an(2,:).*sintheta) + Jc*(bn(3,:).*costheta - an(3,:).*sintheta);

    xN = Geom(par.srcgeom).boundary(1).Xn;
    xL = SimGeom.PlotGeom(par.srcgeom).domain.Xn;

    Xm = Geom(par.srcgeom).boundary(1).Xm;
    Ym = Geom(par.srcgeom).boundary(1).Ym;

    dxm = Xm(2:end) - Xm(1:end-1);
    dym = Ym(2:end) - Ym(1:end-1);

    nx = -dym./sqrt(dxm.^2+dym.^2);
        
    NupN = sum([Geom(1:par.srcgeom).Nbnd]);
    NlowN = NupN + 1 - Geom(par.srcgeom).Nbnd;
    NupL = sum([SimGeom.PlotGeom(1:par.srcgeom).Ldmn]);
    NlowL = NupL + 1 - SimGeom.PlotGeom(par.srcgeom).Ldmn;

    uN(NlowN:NupN) = -(an_nn.*k_nn_2)*sin(k_nn*xN) - (bn_nn.*k_nn_2)*cos(k_nn*xN);
    uL(NlowL:NupL) = -(an_nn.*k_nn_2)*sin(k_nn*xL) - (bn_nn.*k_nn_2)*cos(k_nn*xL);
    qN(NlowN:NupN) = nx.*(-(an_nn.*k_nn_1)*cos(k_nn.*xN) + (bn_nn.*k_nn_1)*sin(k_nn*xN));
      
    xplt = 0:0.001:0.12;
    % xplt = xN;
    figure
    jplt = an_nn*sin(k_nn*xplt) + bn_nn*cos(k_nn*xplt);
    plot(xplt,jplt)
end
