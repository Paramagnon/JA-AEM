function [G00,H00,U00,U00T,PlotGeom,Rngs] = returnExportMatrices(Geom, options, GTlow, GTup, dmnrT)
          
    options.dLx = 0;
    options.dLy = 0;

    options.rmin = Geom(1).domain.rmin/2;
    [PlotGeom(1).domain, PlotGeom(1).Ldmn] = insertDomainPoints(Geom(1), options);

    options.dLx = sqrt(options.dlplot);
    options.dLy = sqrt(options.dlplot);
    for ii = 2:length(Geom)-1
        options.rmin = Geom(ii).domain.rmin/2;
        [PlotGeom(ii).domain, PlotGeom(ii).Ldmn] = insertDomainPoints(Geom(ii), options);
    end

    options.dLx = 0;
    options.dLy = 0;
    
    options.rmin = Geom(end).domain.rmin/2;
    [PlotGeom(end+1).domain, PlotGeom(end+1).Ldmn] = insertDomainPoints(Geom(end), options);
    
    % Find total matrix dimensions    
    N = 0;
    L = 0;
    for kk = 2:length(Geom)-1
        N = N + Geom(kk).Nbnd;
        L = L + Geom(kk).Ldmn;
    end  
    % LT = sum([Geom(GTlow:GTup).LdmnT]);
    % Find total plot matrix dimensions
    M = 0;
    for kk = 1:length(PlotGeom)
        M = M + PlotGeom(kk).Ldmn;
    end

    Rngs.Nrng = sum([Geom(1:GTlow-1).Nbnd])+1 : sum([Geom(1:GTup).Nbnd]);
    Rngs.Lrng = sum([Geom(1:GTlow-1).Ldmn])+1 : sum([Geom(1:GTup).Ldmn]);
    Rngs.Mrng = sum([PlotGeom(1:GTlow-1).Ldmn])+1 : sum([PlotGeom(1:GTup).Ldmn]);
    
    GeomOffsetL = 0;
    GeomOffsetM = 0;
    U00 = zeros(M,L);   
    
    for pp = 1:length(PlotGeom)
        % U functions 
        % , indicates derivative w.r.t. x or y
        % ' indicates derivative w.r.t. r
        % U,xx = U'' (r,x)^2 + U' r,xx
        % U,yy = U'' (r,y)^2 + U' r,yy  
        Lpart = Geom(pp).Ldmn;
        Mpart = PlotGeom(pp).Ldmn;
        
        U00p = zeros(Mpart,Lpart);

        % Calculate U, Uxx, Uyy
        XM = PlotGeom(pp).domain.Xn;
        YM = PlotGeom(pp).domain.Yn;
        XL = Geom(pp).domain.Xn;
        YL = Geom(pp).domain.Yn;
        for ii = 1:PlotGeom(pp).Ldmn    
            x = XM(ii); % Internal node (X)
            y = YM(ii); % Internal node (y)               
            for kk = 1:Geom(pp).Ldmn               
                ksi = XL(kk);
                eta = YL(kk);

                r = sqrt((ksi-x).^2+(eta-y).^2);
                if r < eps
                    r = 0;
                end
                r2 = r^2;

                U00p(ii,kk) = 0.25*r2 + (1/9)*r2*r;
            end
        end
        U00(GeomOffsetM+1:GeomOffsetM+Mpart,GeomOffsetL+1:GeomOffsetL+Lpart) = U00p;
    
        GeomOffsetL = GeomOffsetL + Lpart;
        GeomOffsetM = GeomOffsetM + Mpart;
    end
    % 
    % GeomOffsetLT = 0;
    % GeomOffsetM = 0;
    % U00T = zeros(sum([PlotGeom(GTlow:GTup).Ldmn]),LT);   
    % 
    % for pp = GTlow:GTup
    %     % U functions 
    %     % , indicates derivative w.r.t. x or y
    %     % ' indicates derivative w.r.t. r
    %     % U,xx = U'' (r,x)^2 + U' r,xx
    %     % U,yy = U'' (r,y)^2 + U' r,yy  
    %     LTpart = Geom(pp).LdmnT;
    %     Mpart = PlotGeom(pp).Ldmn;
    % 
    %     U00p = zeros(Mpart,LTpart);
    %     epsilon = dmnrT(pp-GTlow+1);
    %     % Calculate U, Uxx, Uyy
    %     XM = PlotGeom(pp).domain.Xn;
    %     YM = PlotGeom(pp).domain.Yn;
    %     XL = Geom(pp).domainT.Xn;
    %     YL = Geom(pp).domainT.Yn;
    %     for ii = 1:PlotGeom(pp).Ldmn    
    %         x = XM(ii); % Internal node (X)
    %         y = YM(ii); % Internal node (y)               
    %         for kk = 1:Geom(pp).LdmnT
    %             ksi = XL(kk);
    %             eta = YL(kk);
    % 
    %             r = sqrt((ksi-x).^2+(eta-y).^2);
    %             if r < 1e-12
    %                 U00p(ii,kk) = 0;                                               
    %             else
    %                 reps = r/epsilon;
    %                 r2 = reps^2;
    %                 U00p(ii,kk) = 0.25*r2 + (1/9)*r2*reps;
    %             end
    %         end
    %     end
    %     U00T(GeomOffsetM+1:GeomOffsetM+Mpart,GeomOffsetLT+1:GeomOffsetLT+LTpart) = U00p;
    % 
    %     GeomOffsetLT = GeomOffsetLT + LTpart;
    %     GeomOffsetM = GeomOffsetM + Mpart;
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
    GeomOffsetM = 0;

    % Calculate G00, Gxx, Gyy
    G00 = zeros(M,N);    
    % Calculate H00, Hxx, Hyy
    H00 = zeros(M,N);    

    for pp = 1:length(PlotGeom)
        Mpart = PlotGeom(pp).Ldmn;
        Npart = Geom(pp).Nbnd;

        % Calculate G00, Gxx, Gyy
        G00p = zeros(Mpart,Npart);       
        % Calculate H00, Hxx, Hyy
        H00p = zeros(Mpart,Npart);

        XM = PlotGeom(pp).domain.Xn;
        YM = PlotGeom(pp).domain.Yn;

        Noffset = 0;
        for ll=1:length(Geom(pp).boundary)
            Xm = Geom(pp).boundary(ll).Xm; 
            Ym = Geom(pp).boundary(ll).Ym;           
            Nbnd = length(Xm)-1;
    
            G00t = zeros(Mpart,Nbnd);            
            H00t = zeros(Mpart,Nbnd);            
    
            for kk = 1:Mpart
                G00kk = zeros(1,Nbnd);                
                H00kk = zeros(1,Nbnd);
                                
                x = XM(kk);
                y = YM(kk);
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
            
                    if rj < eps
                        G00kk(jj) = (lj/(2*pi))*(log(lj/2)-1);
                        H00kk(jj) = 0;
                    else                       
                        G00kk(jj) = (lj/(4*pi))*(log(rj)*wn);
                        H00kk(jj) = DA/(2*pi);
                    end

                end
                G00t(kk,:) = G00kk; % G00t (Lpart X Nbnd)
                H00t(kk,:) = H00kk; % H00t (Lpart X Nbnd)
            end       
            G00p(:,Noffset+1:Noffset+Nbnd) = G00t;            
            H00p(:,Noffset+1:Noffset+Nbnd) = H00t;

            Noffset = Noffset + Nbnd;
        end

        G00(GeomOffsetM+1:GeomOffsetM+Mpart,GeomOffsetN+1:GeomOffsetN+Npart) = G00p;
        H00(GeomOffsetM+1:GeomOffsetM+Mpart,GeomOffsetN+1:GeomOffsetN+Npart) = H00p;
        
        GeomOffsetM = GeomOffsetM + Mpart;
        GeomOffsetN = GeomOffsetN + Npart;
    end

    if nnz(G00)/numel(G00) < 0.1   
        G00 = sparse(G00);
        H00 = sparse(H00);
    end
    if nnz(U00)/numel(U00) < 0.1   
        U00 = sparse(U00);        
    end
    % if nnz(U00T)/numel(U00T) < 0.1   
    %     U00T = sparse(U00T);        
    % end
    U00T = 0;
end