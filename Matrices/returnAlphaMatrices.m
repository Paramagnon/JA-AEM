function [G00,Gx0,Gy0,H00,Hx0,Hy0,U00,Ux0,Uy0,UxxUyy,Phi] = returnAlphaMatrices(Geom)
    % [G00,Gx0,Gy0,Gxx,Gyy,H00,Hx0,Hy0,Hxx,Hyy,U00,Ux0,Uy0,Uxx,Uyy] 
    
    % LxL matrices for solving alpha each time step, these matrices capture
    % the geometry of the domain so they do not change each time step and
    % are fixed
      
    % Approximation function used: phi = 1 + r 
    % Corresponding particular solution: u = 1/4 r^2 + 1/9 r^3
    
    % Find total matrix dimensions    
    N = 0;
    L = 0;
    for kk = 1:length(Geom)
        N = N + Geom(kk).Nbnd;
        L = L + Geom(kk).Ldmn;
    end  

    GeomOffsetL = 0;

    U00 = zeros(L,L);
    Ux0 = zeros(L,L);
    Uy0 = zeros(L,L);
    Uxx = zeros(L,L);
    Uyy = zeros(L,L);       
    Phi = zeros(L,L);

    for pp = 1:length(Geom)
        % U functions 
        % , indicates derivative w.r.t. x or y
        % ' indicates derivative w.r.t. r
        % U,xx = U'' (r,x)^2 + U' r,xx
        % U,yy = U'' (r,y)^2 + U' r,yy  
        Lpart = Geom(pp).Ldmn;

        U00p = zeros(Lpart,Lpart);
        Ux0p = zeros(Lpart,Lpart);
        Uy0p = zeros(Lpart,Lpart);
        Uxxp = zeros(Lpart,Lpart);
        Uyyp = zeros(Lpart,Lpart);   
        Phip = zeros(Lpart,Lpart);   
        % Calculate U, Uxx, Uyy
        XL = Geom(pp).domain.Xn;
        YL = Geom(pp).domain.Yn;
        for ii = 1:Geom(pp).Ldmn    
            x = XL(ii); % Internal node (X)
            y = YL(ii); % Internal node (y)               
            for kk = 1:Geom(pp).Ldmn
                if kk ~= ii
                    ksi = XL(kk);
                    eta = YL(kk);

                    r = sqrt((ksi-x).^2+(eta-y).^2);
                    if r < eps
                        U00p(ii,kk) = 0;
                        Ux0p(ii,kk) = 0;
                        Uy0p(ii,kk) = 0;
                        Uxxp(ii,kk) = 0.5;
                        Uyyp(ii,kk) = 0.5;
                        Phip(ii,kk) = 1;
                    else
                        r2 = r^2;

                        rx = (ksi-x)/r;             
                        ry = (eta-y)/r;

                        rxx = (ry^2)/r;
                        ryy = (rx^2)/r;

                        U_ = (0.5*r + (1/3)*r2);
                        U__ = (0.5 + (2/3)*r);

                        U00p(ii,kk) = 0.25*r2 + (1/9)*r2*r;                    
                        Ux0p(ii,kk) = rx*(0.5*r + (1/3)*r2);
                        Uy0p(ii,kk) = ry*(0.5*r + (1/3)*r2);                    
                        Uxxp(ii,kk) = U__* rx^2 + U_*rxx;
                        Uyyp(ii,kk) = U__* ry^2 + U_*ryy;
                        Phip(ii,kk) = 1 + r;
                    end

                else
                    U00p(ii,kk) = 0;
                    Ux0p(ii,kk) = 0;
                    Uy0p(ii,kk) = 0;
                    Uxxp(ii,kk) = 0.5;
                    Uyyp(ii,kk) = 0.5;
                    Phip(ii,kk) = 1;
                end           
            end
        end
        U00(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = U00p;
        Ux0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Ux0p;
        Uy0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Uy0p;
        Uxx(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Uxxp;
        Uyy(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Uyyp;
        Phi(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Phip;

        GeomOffsetL = GeomOffsetL + Lpart;
    end

    % epsilon = SimGeom.epsilon;     
    % for pp = 1:length(Geom)
    %     Lpart = Geom(pp).Ldmn;
    % 
    %     U00p = zeros(Lpart,Lpart);
    %     Ux0p = zeros(Lpart,Lpart);
    %     Uy0p = zeros(Lpart,Lpart);
    %     Uxxp = zeros(Lpart,Lpart);
    %     Uyyp = zeros(Lpart,Lpart);   
    %     Phip = zeros(Lpart,Lpart);   
    % 
    %     XL = Geom(pp).domain.Xn;
    %     YL = Geom(pp).domain.Yn;
    % 
    %     for ii = 1:Lpart    
    %         x = XL(ii);
    %         y = YL(ii);
    % 
    %         for kk = 1:Lpart
    %             if kk ~= ii
    %                 ksi = XL(kk);
    %                 eta = YL(kk);
    % 
    %                 dx = ksi - x;
    %                 dy = eta - y;
    %                 r2 = dx^2 + dy^2;
    %                 r  = sqrt(r2);
    % 
    %                 if r < 1e-12
    %                     U00p(ii,kk) = 0;
    %                     Ux0p(ii,kk) = 0;
    %                     Uy0p(ii,kk) = 0;
    %                     Uxxp(ii,kk) = 0;
    %                     Uyyp(ii,kk) = 0;
    %                     Phip(ii,kk) = 1;
    %                 else
    %                     er = exp(-(epsilon*r)^2);
    %                     fac1 = -0.5 * er / r2;
    %                     fac2 = er / (2*r2^2);
    % 
    %                     U00p(ii,kk) = -1/(4*epsilon^2) * expint(-(epsilon*r)^2);
    %                     Ux0p(ii,kk) = fac1 * dx;
    %                     Uy0p(ii,kk) = fac1 * dy;
    %                     Uxxp(ii,kk) = fac2 * (2*dx^2 - r2);
    %                     Uyyp(ii,kk) = fac2 * (2*dy^2 - r2);
    %                     Phip(ii,kk) = er;
    %                 end
    % 
    %             else                                        
    %                 U00p(ii,kk) = 0;
    %                 Ux0p(ii,kk) = 0;
    %                 Uy0p(ii,kk) = 0;                   
    %                 Uxxp(ii,kk) = -epsilon^2/2;
    %                 Uyyp(ii,kk) = -epsilon^2/2;
    %                 Phip(ii,kk) = 1;
    %             end
    %         end
    %     end
    %     U00(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = U00p;
    %     Ux0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Ux0p;
    %     Uy0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Uy0p;
    %     Uxx(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Uxxp;
    %     Uyy(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Uyyp;
    %     Phi(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetL+1:GeomOffsetL+Lpart) = Phip;
    % 
    %     GeomOffsetL = GeomOffsetL + Lpart;
    % end

    % Gaussian quadrature 
    %--- n = 4 coefficients
    % ksi = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053]; % Gaussian quadrature integration points, from table for n = 4
    % wn = [0.347854845137454; 0.652145154862546; 0.652145154862546; 0.347854845137454]; % Gaussian quadrature integration weights, from table for n = 4

    % Gaussian quadrature
    %--- n = 7 coefficients
    % ksi = [-0.949107912342759, -0.741531185599394, -0.405845151377397, 0.000000000000000, 0.405845151377397, 0.741531185599394, 0.949107912342759];
    % wn = [0.129484966168870; 0.279705391489277; 0.381830050505119; 0.417959183673469; 0.381830050505119; 0.279705391489277; 0.129484966168870];

    % Gaussian quadrature
    %--- n = 16 coefficients
    % ksi = [-0.989400934991649, -0.944575023073232, -0.865631202387831, -0.755404408355003, -0.617876244402643, -0.458016777657227, -0.281603550779258, -0.095012509837637, 0.095012509837637, 0.281603550779258, 0.458016777657227, 0.617876244402643, 0.755404408355003, 0.865631202387831, 0.944575023073232, 0.989400934991649];
    % wn = [0.027152459411754; 0.062253523938647; 0.095158511682492; 0.124628971255533; 0.149595988816576; 0.169156519395002; 0.182603415044923; 0.189450610455068; 0.189450610455068; 0.182603415044923; 0.169156519395002; 0.149595988816576; 0.124628971255533; 0.095158511682492; 0.062253523938647; 0.027152459411754];
    
    % Gaussian quadrature
    % --- n = 36 coefficients
    % ksi = [-0.0430181984737086, 0.0430181984737086, -0.1287361038093848, 0.1287361038093848, -0.2135008923168656, 0.2135008923168656, -0.2966849953440283, 0.2966849953440283, -0.3776725471196892, 0.3776725471196892,-0.4558639444334203, 0.4558639444334203,-0.5306802859262452, 0.5306802859262452, -0.6015676581359806, 0.6015676581359806, -0.6680012365855210, 0.6680012365855210, -0.7294891715935565, 0.7294891715935565, -0.7855762301322066, 0.7855762301322066, -0.8358471669924753, 0.8358471669924753, -0.8799298008903972, 0.8799298008903972, -0.9174977745156591, 0.9174977745156591, -0.9482729843995076, 0.9482729843995076, -0.9720276910496980, 0.9720276910496980, -0.9885864789022122, 0.9885864789022122, -0.9978304624840858, 0.9978304624840858]; % Gaussian quadrature integration points, from table for n = 36
    % wn = [0.0859832756703948; 0.0859832756703948; 0.0853466857393386; 0.0853466857393386; 0.0840782189796619; 0.0840782189796619; 0.0821872667043397; 0.0821872667043397; 0.0796878289120716; 0.0796878289120716; 0.0765984106458707; 0.0765984106458707; 0.0729418850056531; 0.0729418850056531; 0.0687453238357364; 0.0687453238357364; 0.0640397973550155; 0.0640397973550155; 0.0588601442453248; 0.0588601442453248; 0.0532447139777599; 0.0532447139777599; 0.0472350834902660; 0.0472350834902660; 0.0408757509236449; 0.0408757509236449; 0.0342138107703072; 0.0342138107703072; 0.0272986214985688; 0.0272986214985688; 0.0201815152977355; 0.0201815152977355; 0.0129159472840656; 0.0129159472840656; 0.0055657196642450; 0.0055657196642450]; % Gaussian quadrature integration weights, from table for n = 36;

    % Gaussian quadrature
    % --- n = 64 coefficients
    ksi = [-0.0243502926634244,0.0243502926634244,-0.0729931217877990,0.0729931217877990,-0.1214628192961206,0.1214628192961206,-0.1696444204239928,0.1696444204239928,-0.2174236437400071,0.2174236437400071,-0.2646871622087674,0.2646871622087674,-0.3113228719902110,0.3113228719902110,-0.3572201583376681,0.3572201583376681,-0.4022701579639916,0.4022701579639916,-0.4463660172534641,0.4463660172534641,-0.4894031457070530,0.4894031457070530,-0.5312794640198946,0.5312794640198946,-0.5718956462026340,0.5718956462026340,-0.6111553551723933,0.6111553551723933,-0.6489654712546573,0.6489654712546573,-0.6852363130542333,0.6852363130542333,-0.7198818501716109,0.7198818501716109,-0.7528199072605319,0.7528199072605319,-0.7839723589433414, ...
        0.7839723589433414,-0.8132653151227975,0.8132653151227975,-0.8406292962525803,0.8406292962525803,-0.8659993981540928,0.8659993981540928,-0.8893154459951141,0.8893154459951141,-0.9105221370785028,0.9105221370785028,-0.9295691721319396,0.9295691721319396,-0.9464113748584028,0.9464113748584028,-0.9610087996520538,0.9610087996520538,-0.9733268277899110,0.9733268277899110,-0.9833362538846260,0.9833362538846260,-0.9910133714767443,0.9910133714767443,-0.9963401167719553,0.9963401167719553,-0.9993050417357722,0.9993050417357722]; % Gaussian quadrature integration points, from table for n = 64 
    wn = [0.0486909570091397;0.0486909570091397;0.0485754674415034;0.0485754674415034;0.0483447622348030;0.0483447622348030;0.0479993885964583;0.0479993885964583;0.0475401657148303;0.0475401657148303;0.0469681828162100;0.0469681828162100;0.0462847965813144;0.0462847965813144;0.0454916279274181;0.0454916279274181;0.0445905581637566;0.0445905581637566;0.0435837245293235;0.0435837245293235;0.0424735151236536;0.0424735151236536;0.0412625632426235;0.0412625632426235;0.0399537411327203;0.0399537411327203;0.0385501531786156;0.0385501531786156;0.0370551285402400;0.0370551285402400;0.0354722132568824;0.0354722132568824;0.0338051618371416;0.0338051618371416;0.0320579283548516;0.0320579283548516;0.0302346570724025;0.0302346570724025;...
        0.0283396726142595;0.0283396726142595;0.0263774697150547;0.0263774697150547;0.0243527025687109;0.0243527025687109;0.0222701738083833;0.0222701738083833;0.0201348231535302;0.0201348231535302;0.0179517157756973;0.0179517157756973;0.0157260304760247;0.0157260304760247;0.0134630478967186;0.0134630478967186;0.0111681394601311;0.0111681394601311;0.0088467598263639;0.0088467598263639;0.0065044579689784;0.0065044579689784;0.0041470332605625;0.0041470332605625;0.0017832807216964;0.0017832807216964]; % Gaussian quadrature integration weights, from table for n = 64

    GeomOffsetN = 0;
    GeomOffsetL = 0;

    % Calculate G00, Gxx, Gyy
    G00 = zeros(L,N);
    Gx0 = zeros(L,N);
    Gy0 = zeros(L,N);
    % Gxx = zeros(L,N);
    % Gyy = zeros(L,N);
    % Calculate H00, Hxx, Hyy
    H00 = zeros(L,N);
    Hx0 = zeros(L,N);
    Hy0 = zeros(L,N);
    % Hxx = zeros(L,N);
    % Hyy = zeros(L,N);

    for pp = 1:length(Geom)
        Lpart = Geom(pp).Ldmn;
        Npart = Geom(pp).Nbnd;

        % Calculate G00, Gxx, Gyy
        G00p = zeros(Lpart,Npart);
        Gx0p = zeros(Lpart,Npart);
        Gy0p = zeros(Lpart,Npart);
        % Gxxp = zeros(Lpart,Npart);
        % Gyyp = zeros(Lpart,Npart);

        % Calculate H00, Hxx, Hyy
        H00p = zeros(Lpart,Npart);
        Hx0p = zeros(Lpart,Npart);
        Hy0p = zeros(Lpart,Npart);
        % Hxxp = zeros(Lpart,Npart);
        % Hyyp = zeros(Lpart,Npart);

        XL = Geom(pp).domain.Xn;
        YL = Geom(pp).domain.Yn;

        Noffset = 0;
        for ll=1:length(Geom(pp).boundary)
            Xm = Geom(pp).boundary(ll).Xm; 
            Ym = Geom(pp).boundary(ll).Ym;           
            Nbnd = length(Xm)-1;
    
            G00t = zeros(Lpart,Nbnd);
            Gx0t = zeros(Lpart,Nbnd);
            Gy0t = zeros(Lpart,Nbnd);
            % Gxxt = zeros(Lpart,Nbnd);
            % Gyyt = zeros(Lpart,Nbnd);
            H00t = zeros(Lpart,Nbnd);
            Hx0t = zeros(Lpart,Nbnd);
            Hy0t = zeros(Lpart,Nbnd);
            % Hxxt = zeros(Lpart,Nbnd);
            % Hyyt = zeros(Lpart,Nbnd);
    
            for kk = 1:Lpart
                G00kk = zeros(1,Nbnd);
                Gx0kk = zeros(1,Nbnd);
                Gy0kk = zeros(1,Nbnd);
                % Gxxkk = zeros(1,Nbnd);
                % Gyykk = zeros(1,Nbnd);
                H00kk = zeros(1,Nbnd);
                Hx0kk = zeros(1,Nbnd);
                Hy0kk = zeros(1,Nbnd);
                % Hxxkk = zeros(1,Nbnd);
                % Hyykk = zeros(1,Nbnd);
                
                x = XL(kk);
                y = YL(kk);
                for jj = 1:Nbnd
                    lj = sqrt((Xm(jj+1)-Xm(jj))^2+(Ym(jj+1)-Ym(jj))^2); % Length of elements, use extreme values
                    
                    dy1 = (Ym(jj)-y);
                    dx1 = (Xm(jj)-x);
                    dy2 = (Ym(jj+1)-y);
                    dx2 = (Xm(jj+1)-x);
    
                    DL1 = sqrt(dx1^2+dy1^2);
                    COS1 = dx1/DL1;
                    SIN1 = dy1/DL1;
        
                    DX2R = dx2*COS1 + dy2*SIN1;
                    DY2R = -dx2*SIN1 + dy2*COS1;
                    DA = atan2(DY2R,DX2R);
        
                    x_ksi0 = (Xm(jj+1)+Xm(jj))/2;
                    x_ksi1 = (Xm(jj+1)-Xm(jj))/2;
                    y_ksi0 = (Ym(jj+1)+Ym(jj))/2;
                    y_ksi1 = (Ym(jj+1)-Ym(jj))/2;
    
                    x_ksi = x_ksi0 + x_ksi1*ksi;
                    y_ksi = y_ksi0 + y_ksi1*ksi;
    
                    rj = sqrt((x_ksi-x).^2 + (y_ksi-y).^2);
    
                    drx2 = (x_ksi - x).^2; 
                    dry2 = (y_ksi - y).^2;
                   
                    drx = (x_ksi - x);
                    dry = (y_ksi - y);
                    
                    rj2p = rj.^2;

                    rj2 = 1./(2*pi*rj2p);
                    rj4 = rj2./rj2p;
                    % rj6 = rj4./rj2p;
    
                    dxm = Xm(jj+1) - Xm(jj);
                    dym = Ym(jj+1) - Ym(jj);
    
                    % nx = -dxm/sqrt(dxm^2+dym^2);
                    % ny = dym/sqrt(dxm^2+dym^2);
                    nx = dym/sqrt(dxm^2+dym^2);
                    ny = -dxm/sqrt(dxm^2+dym^2);

                    % Hfac = -5*dry.*drx2*ny + 2*dry.*dry2*ny + ...
                    %     5*drx.*dry2*nx - 2*drx.*drx2*nx;
                    
                    lj_2 = lj/2;
                        
                    %--- Experimental
                    if (min(rj)>1e-12)
                        G00kk(jj) = (lj/(4*pi))*(log(rj)*wn);
                        Gx0kk(jj) = lj_2*(rj2.*drx)*wn;
                        Gy0kk(jj) = lj_2*(rj2.*dry)*wn;
                        % Gxxkk(jj) = lj_2*(rj4.*(dry2-drx2))*wn;
                        % Gyykk(jj) = -Gxxkk(jj); % (rj4.*(drx2-dry2))*wn
        
                        H00kk(jj) = DA/(2*pi);
                        Hx0kk(jj) = lj_2*(rj4.*((drx2-dry2)*nx+2.*drx.*dry*ny))*wn;
                        Hy0kk(jj) = lj_2*(rj4.*((dry2-drx2)*ny+2.*drx.*dry*nx))*wn;
                        % Hxxkk(jj) = lj_2*(rj6.*(-Hfac)*wn);
                        % Hyykk(jj) = -Hxxkk(jj); % (rj6.*Hfac*wn)
                    else
                        G00kk(jj) = (lj/(2*pi))*(log(lj/2)-1);
                        Gx0kk(jj) = 0;
                        Gy0kk(jj) = 0;
                        H00kk(jj) = -0.5;
                        Hx0kk(jj) = 0;
                        Hy0kk(jj) = 0;
                    end

                end
                G00t(kk,:) = G00kk; % G00t (Lpart X Nbnd)
                Gx0t(kk,:) = Gx0kk; % Gx0t (Lpart X Nbnd)
                Gy0t(kk,:) = Gy0kk; % Gy0t (Lpart X Nbnd)
                % Gxxt(kk,:) = Gxxkk; % Gxxt (Lpart X Nbnd) 
                % Gyyt(kk,:) = Gyykk; % Gyyt (Lpart X Nbnd)
                H00t(kk,:) = H00kk; % H00t (Lpart X Nbnd)
                Hx0t(kk,:) = Hx0kk; % Hx0t (Lpart X Nbnd)
                Hy0t(kk,:) = Hy0kk; % Hy0t (Lpart X Nbnd)
                % Hxxt(kk,:) = Hxxkk; % Hxxt (Lpart X Nbnd)
                % Hyyt(kk,:) = Hyykk; % Hyyt (Lpart X Nbnd)
            end       
            G00p(:,Noffset+1:Noffset+Nbnd) = G00t;
            Gx0p(:,Noffset+1:Noffset+Nbnd) = Gx0t;
            Gy0p(:,Noffset+1:Noffset+Nbnd) = Gy0t;
            % Gxxp(:,Noffset+1:Noffset+Nbnd) = Gxxt;
            % Gyyp(:,Noffset+1:Noffset+Nbnd) = Gyyt;
            H00p(:,Noffset+1:Noffset+Nbnd) = H00t;
            Hx0p(:,Noffset+1:Noffset+Nbnd) = Hx0t;
            Hy0p(:,Noffset+1:Noffset+Nbnd) = Hy0t;
            % Hxxp(:,Noffset+1:Noffset+Nbnd) = Hxxt;
            % Hyyp(:,Noffset+1:Noffset+Nbnd) = Hyyt;
            Noffset = Noffset + Nbnd;
        end

        G00(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = G00p;
        Gx0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Gx0p;
        Gy0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Gy0p;
        % Gxx(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Gxxp;
        % Gyy(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Gyyp;
        H00(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = H00p;
        Hx0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Hx0p;
        Hy0(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Hy0p;
        % Hxx(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Hxxp;
        % Hyy(GeomOffsetL+1:GeomOffsetL+Lpart,GeomOffsetN+1:GeomOffsetN+Npart) = Hyyp;

        GeomOffsetL = GeomOffsetL + Lpart;
        GeomOffsetN = GeomOffsetN + Npart;
    end

    if nnz(G00)/numel(G00) < 0.1
        G00 = sparse(G00);
        Gx0 = sparse(Gx0);
        Gy0 = sparse(Gy0);
        H00 = sparse(H00);
        Hx0 = sparse(Hx0);
        Hy0 = sparse(Hy0);
    end
    UxxUyy = Uxx+Uyy;
end