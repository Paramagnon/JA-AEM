function [Geom] = linkNghbrBndrs(Geom) 
    
    for kk = 1:length(Geom)
        N = 0;
        Nbnd = zeros(1,length(Geom(kk).boundary));
        for nn = 1:length(Geom(kk).boundary)
            Nbnd(nn) = sum(Geom(kk).boundary(nn).Nidx);
            N = N + Nbnd(nn);
        end
        Geom(kk).idx_linked = zeros(2,N);   % Store linking idx [idx; domain]
        Geom(kk).Nbnd = Nbnd;               % Length of each boundary inside geom
        Geom(kk).SC = [];
        % Geom(kk).u_sol = zeros(1,N);
        % Geom(kk).un_sol = zeros(1,N);       
        % Geom(kk).u_p = zeros(1,N);          % Particular solution for domain source
        % Geom(kk).un_p = zeros(1,N);         % Particular solution for domain source
    end
       
    for kk = 1:length(Geom)
        Noffset = 0;
        SCit = 1;
        for nn = 1:length(Geom(kk).boundary)            
            for ll = 1:length(Geom(kk).boundary(nn).Nidx)               
                domain_idx = Geom(kk).boundary(nn).domain_idx(ll);
                Nbnd = Geom(kk).boundary(nn).Nidx(ll);
                if domain_idx ~= 0
                    boundary_idx = Geom(kk).boundary(nn).boundary_idx(ll);
                    element_idx = Geom(kk).boundary(nn).element_idx(ll);
    
                    % Fetch indices in other geom/ boundary  
                    Nbnd_offset = sum(Geom(domain_idx).Nbnd(1:boundary_idx-1));
                    Nstart = Nbnd_offset + sum(Geom(domain_idx).boundary(boundary_idx).Nidx(1:element_idx-1)) + 1;
                    Nend = Nstart + Geom(domain_idx).boundary(boundary_idx).Nidx(element_idx) - 1;
                    
                    % Fetch length of boundary                   
                    Geom(kk).idx_linked(1,1+Noffset:Noffset+Nbnd) = Nend:-1:Nstart;
                    Geom(kk).idx_linked(2,1+Noffset:Noffset+Nbnd) = domain_idx; 

                    Nbndr_idx = sum(Geom(kk).boundary(nn).Nidx(1:ll-1))+1;
                    if Geom(kk).boundary(nn).bndr_idx(Nbndr_idx) == 4
                        Geom(kk).SC(end+1).idx = 1+Noffset:Noffset+Nbnd;
                        Nscrange = Nbndr_idx+1:Nbndr_idx+Nbnd;
                        Geom(kk).SC(end).dy = abs(Geom(kk).boundary(nn).Ym(Nscrange)-Geom(kk).boundary(nn).Ym(Nscrange-1));
                        Geom(kk).SC(end).pol = Geom(kk).SCpol(SCit);
                        SCit = SCit + 1;
                    end
                end
                Noffset = Noffset + Nbnd;
            end
        end        
    end

    for kk = 1:length(Geom)
        Geom(kk).Nbnd = sum(Geom(kk).Nbnd);
    end
end




