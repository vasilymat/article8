function [ f ] = plh( dz, kz, r, T, lam, dx )

if (lam*dz/(size(T,1)*dx^2)) < 1 
    f = iFTS( FTS(T).*(exp(1j*kz.*dz)) );
else
    k = 2*pi/lam;
    fh = exp(1j*k.*dz)./(1j*lam.*dz)*exp(1j*k.*(r.^2)/2./dz);
    fh = fh.*dx^2;
    f = iFTS( FTS(T).*FTS(fh) );
end




    