%%
L = Geom(1).Ldmn + Geom(2).Ldmn;
phi = zeros(L,L);
Loffset = 0;

for pp = 1:length(Geom)
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
                phi(Loffset+ii,Loffset+kk) = 1+r;
            end
        end
    end
    Loffset = Loffset + Geom(pp).Ldmn;
end