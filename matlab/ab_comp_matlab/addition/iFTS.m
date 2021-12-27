function [F] = iFTS(f)
%shift anf Fourier transform
%   Detailed explanation goes here
F = ifftshift(ifft2(ifftshift(f)));
end

