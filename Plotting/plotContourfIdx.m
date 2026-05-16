function [idx, xknwn, yknwn, xp, yp] = plotContourfIdx(Geom, dlplot)
    % rmin = 1e-6; % Set small since only for plotting

    xknwn = zeros(1,Geom.Nbnd+Geom.Ldmn);
    yknwn = zeros(1,Geom.Nbnd+Geom.Ldmn);    

    Noffset = 0;
    for ii = 1:length(Geom.boundary)
        Npart = length(Geom.boundary(ii).bndr_idx);

        xknwn(Noffset+1:Noffset+Npart) = Geom.boundary(ii).Xn;
        yknwn(Noffset+1:Noffset+Npart) = Geom.boundary(ii).Yn;
        Noffset = Noffset + Npart;
    end
    xknwn(Noffset+1:Noffset+Geom.Ldmn) = Geom.domain.Xn;
    yknwn(Noffset+1:Noffset+Geom.Ldmn) = Geom.domain.Yn;

    xmin = min(xknwn); 
    xmax = max(xknwn);
    ymin = min(yknwn); 
    ymax = max(yknwn);

    xi = linspace(xmin,xmax,round((xmax-xmin)*sqrt(dlplot),0));
    yi = linspace(ymin,ymax,round((ymax-ymin)*sqrt(dlplot),0));
    
    [Xi,Yi] = meshgrid(xi, yi);   

    % Find nodes that are outside the boundaries, use the idea of ray
    % casting and seeing through how many edges the ray travels, if it is
    % uneven the point inside and vice verse. Draw line parallel to x-axis
    % from y point of the node, simplifies checks 

    idx1 = ones(size(Xi));
    Boundary = Geom.boundary;

    parfor ii = 2:size(Yi,1)-1
        idx0 = ones(1,size(Yi,2));
        for ll = 2:size(Yi,2)-1
            passes = 0;
            for kk = 1:length(Boundary)
                for nn = 1:length(Boundary(kk).Ym)-1
                    Y1 = Boundary(kk).Ym(nn);
                    Y2 = Boundary(kk).Ym(nn+1);                
                    dY = Y2-Y1;
                    if Yi(ii,ll) > min([Y1,Y2]) && Yi(ii,ll) <= max([Y1,Y2]) && dY ~= 0
                        X1 = Boundary(kk).Xm(nn);
                        X2 = Boundary(kk).Xm(nn+1);
                        Xint = (Yi(ii,ll)-Y1)*(X2-X1)/(Y2-Y1) + X1;
                        if Xi(ii,ll) > Xint
                            passes = passes+1;
                        end
                    end
                end
            end
            idx0(ll) = mod(passes,2); 
        end
        idx1(ii,:) = idx0; 
    end
    
    % Find nodes which do not satisfy the min distance condition to the
    % boundary       
    % idx2 = zeros(size(Xi));
    % parfor ii = 1:size(Yi,1)
    %     idx0 = zeros(1,size(Yi,2));
    %     for ll = 1:size(Yi,2)
    %         for kk = 1:length(Boundary)
    %             for nn = 1:length(Boundary(kk).Ym)-1
    %                 Y1 = Boundary(kk).Ym(nn);
    %                 Y2 = Boundary(kk).Ym(nn+1);
    %                 X1 = Boundary(kk).Xm(nn);
    %                 X2 = Boundary(kk).Xm(nn+1);
    % 
    %                 A = Xi(ii,ll) - X1;
    %                 B = Yi(ii,ll) - Y1;
    %                 C = X2 - X1;
    %                 D = Y2 - Y1;
    % 
    %                 dist_case = (A*C + B*D)/(C*C + D*D);
    %                 if dist_case < 0
    %                     xx = X1;
    %                     yy = Y1;
    %                 elseif dist_case > 1
    %                     xx = X2;
    %                     yy = Y2;
    %                 else
    %                     xx = X1 + dist_case*C;
    %                     yy = Y1 + dist_case*D;
    %                 end
    % 
    %                 dx = Xi(ii,ll) - xx;
    %                 dy = Yi(ii,ll) - yy;
    %                 r = sqrt(dx*dx + dy*dy);
    %                 if r >= rmin
    %                     idx0(ll) = 1;                
    %                 end           
    %             end
    %         end
    %     end
    %     idx2(ii,:) = idx0;
    % end
    idx = ~idx1;
    xp = Xi;
    yp = Yi;
end