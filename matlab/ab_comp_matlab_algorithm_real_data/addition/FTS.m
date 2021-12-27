function [F] = FTS(f)
%shift anf Fourier transform
%   Detailed explanation goes here
F = fftshift(fft2(fftshift(f)));
end

