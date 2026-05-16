function f = currentSine(Geom,par,t)
    Jref = par.Jref;
    F = par.F;
    
    f = zeros(sum([Geom.Ldmn]),1);
    % f(Geom(1).Ldmn+1:end) = -Geom(2).mu_f(270)*ones(Geom(2).Ldmn,1)*Jref*sin(2*pi*F*t);
    f(Geom(1).Ldmn+1:end) = ones(Geom(2).Ldmn,1)*Jref*sin(2*pi*F*t);

end
