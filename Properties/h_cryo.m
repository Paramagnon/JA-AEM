function h = h_cryo(options,par)
    A = pi*par.r_pipes^2;                   % [m^2]
    u = 0.001*options.cryoflow/A;           % [m/s]
    L = 2*par.r_pipes;                      % [m] Diameter of the pipes
    Tcryo = par.Tcryo; 
    
    switch options.cryotype
        case 1 % Liquid Nitrogen
            if Tcryo > 77.355 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] above boiling point'])
            elseif Tcryo < 63.151 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] below data range'])
            end
            Pr = Cp_lnitrogen(Tcryo)*muv_lnitrogen(Tcryo)/kappa_lnitrogen(Tcryo);
            Re = rho_lnitrogen(Tcryo)*u*L/muv_lnitrogen(Tcryo);
            Nu = 0.023*Re^0.8*Pr^0.4; % --> Check experimental validation ranges
            h = Nu*kappa_lnitrogen(Tcryo)/L;
  
        case 2 % Liquid Hydrogen           
            if Tcryo > 20.369 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] above boiling point'])
            elseif Tcryo < 13.957 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] below data range'])
            end
            Pr = Cp_lhydrogen(Tcryo)*muv_lhydrogen(Tcryo)/kappa_lhydrogen(Tcryo);
            Re = rho_lhydrogen(Tcryo)*u*L/muv_lhydrogen(Tcryo);
            Nu = 0.023*Re^0.8*Pr^0.4;
            h = Nu*kappa_lhydrogen(Tcryo)/L;

        case 3 % Liquid Helium            
            if Tcryo > 4.2238 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] above boiling point'])
            elseif Tcryo < 2.177 %[K]
                error(['Cryo temp: ' num2str(Tcryo) ' [K] below data range'])
            end
            Pr = Cp_lhelium(Tcryo)*muv_lhelium(Tcryo)/kappa_lhelium(Tcryo);
            Re = rho_lhelium(Tcryo)*u*L/muv_lhelium(Tcryo);
            Nu = 0.022*Re^0.8*Pr^0.4;
            h = Nu*kappa_lhelium(Tcryo)/L;

        otherwise
            disp("Not a valid cryogen type");
    end
end