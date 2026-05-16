function plotMesh(Geom,ifarrow)
    figure
    hold on 
    title("Magnetic domain mesh")
    if isfield(Geom(1),"boundary")
        for kk = 1:length(Geom)
            plot(Geom(kk).domain.Xn,Geom(kk).domain.Yn,".b","MarkerSize",8)
            for nn = 1:length(Geom(kk).boundary)
                idx0 = Geom(kk).boundary(nn).bndr_idx == 0;
                idx1 = Geom(kk).boundary(nn).bndr_idx == 1;
                idx2 = Geom(kk).boundary(nn).bndr_idx == 2;
                idx3 = Geom(kk).boundary(nn).bndr_idx == 3;
                idx4 = Geom(kk).boundary(nn).bndr_idx == 4;
        
                plot(Geom(kk).boundary(nn).Xm(idx0),Geom(kk).boundary(nn).Ym(idx0),"xk","MarkerSize",8)
                plot(Geom(kk).boundary(nn).Xm(idx1),Geom(kk).boundary(nn).Ym(idx1),"xr","MarkerSize",8)
                plot(Geom(kk).boundary(nn).Xm(idx2),Geom(kk).boundary(nn).Ym(idx2),"xg","MarkerSize",8)
                plot(Geom(kk).boundary(nn).Xm(idx3),Geom(kk).boundary(nn).Ym(idx3),"xm","MarkerSize",8)
                plot(Geom(kk).boundary(nn).Xm(idx4),Geom(kk).boundary(nn).Ym(idx4),"xc","MarkerSize",8)
                if ifarrow
                    addArrowLine(Geom(kk).boundary(nn).xgeom,Geom(kk).boundary(nn).ygeom)    
                end
            end
        end
    else
        for kk = 1:length(Geom)
            plot(Geom(kk).domain.Xn,Geom(kk).domain.Yn,".b","MarkerSize",8)
        end
    end
    grid on
    xlabel("x-axis [m]","Interpreter","latex","FontSize",12)
    ylabel("y-axis [m]","Interpreter","latex","FontSize",12)
end