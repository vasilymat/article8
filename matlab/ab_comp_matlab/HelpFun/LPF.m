function [ Ap ] = LPF( Ap, B )
%‘ильтр низких частот. ћножит на гаусс шириной L/B
L = length(Ap);
[x,y] = meshgrid(-L/2:1:L/2-1);
is_in_circle_t = circ(sqrt(x.^2+y.^2)/B); %круг заданной ширины

sig = B;
Pst = exp( -1*( x.^2+y.^2 )/(2*sig^2) )/sqrt(2*pi)/sig ; %ядро свертки дл€ круга

% Psm = ifft2(fft2(fftshift(Pst)).*fft2(is_in_circle_t)); %—ворачиваем и нормируем
% m_psm = max(max(Psm));
% Psm = Psm/m_psm;
% 
% Ap = ifft2( fft2(Ap).*fftshift(Psm) );

Ap = ifft2( fft2(fftshift(Ap)).*fft2(Pst) );
mm = max(max(abs(Ap)));

Ap = abs(Ap)/mm;


end

