function [dmnr,dmnrT] = calculateDomainNormalization(Geom)    
    % Thermal Geom indices
    GTlow = Geom(1).Tlow;
    GTup = Geom(1).Tup;  

    % Normalization per domain
    dmnr = zeros(1,length(Geom));
    dmnrT = zeros(1,GTup-GTlow);
    
    % Distances    
    for pp = 1:length(Geom)  
        r = zeros(1,Geom(pp).Ldmn*Geom(pp).Nbnd);
        idx = 1;
        for ii = 1:Geom(pp).Ldmn
            x = Geom(pp).domain.Xn(ii); % Internal node (y)
            y = Geom(pp).domain.Yn(ii); % Internal node (y)
                        
            for kk = 1:length(Geom(pp).boundary)
                Xn = Geom(pp).boundary(kk).Xn;
                Yn = Geom(pp).boundary(kk).Yn;

                for mm = 1:length(Geom(pp).boundary(kk).Xn)
                    r(idx) = sqrt((x-Xn(mm))^2+(y-Yn(mm))^2);
                    idx = idx + 1;
                end
            end
        end       
        dmnr(pp) = mean(r);
    end

    % Distances   
    for pp = GTlow:GTup
        r = zeros(1,Geom(pp).LdmnT*Geom(pp).Nbnd);
        idx = 1;
        for ii = 1:Geom(pp).LdmnT
            x = Geom(pp).domainT.Xn(ii); % Internal node (y)
            y = Geom(pp).domainT.Yn(ii); % Internal node (y)

            for kk = 1:length(Geom(pp).boundary)
                Xn = Geom(pp).boundary(kk).Xn;
                Yn = Geom(pp).boundary(kk).Yn;

                for mm = 1:length(Geom(pp).boundary(kk).Xn)
                    r(idx) = sqrt((x-Xn(mm))^2+(y-Yn(mm))^2);
                    idx = idx + 1;
                end
            end
        end
        dmnrT(pp-GTlow+1) = mean(r);
    end
end

