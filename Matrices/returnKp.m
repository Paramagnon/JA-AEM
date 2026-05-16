function [K,p] = returnKp(Geom)
    
    % K can be adjusted to include more spatial derivatives
    % K = A Sxx + C Syy + F S00
    % K = Geom.A * Geom.Sxx + Geom.C * Geom.Syy; % Geom.F * Geom.S00; No A dependence
    %--- Simplifies since Hxx = -Hyy and Gxx = -Gyy
    K = Geom.Uxx + Geom.Uyy; % Geom.F * Geom.S00; No A dependence

    % p should be adjusted if the boundary support has time derivatives
    % p = f - (A zxx + C zyy + F z00)
    % p = Geom.f - (Geom.A * Geom.zxx + Geom.C * Geom.zyy); % + Geom.sigma * Geom.dAdt; % No A dependence    
    p = Geom.f; % + Geom.sigma * Geom.dAdt; % No A dependence    
end