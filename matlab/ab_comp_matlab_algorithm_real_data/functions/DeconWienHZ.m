function [Dc,Phase] = DeconWienHZ(ATF,Img,noise_area,Deg)
% Decomposition amplitude ATF by Gauss-Hermite pol. and phase by Zernike
% pol.
% ATF - 64x64 ATF
% Img - aberrated image
% noise_area - is area by with noise is calculated
% H64 and H256 - are arrays of G-H polynomials on the correspondings grids

dx = 300/256*10^-3;
C = size(Img,1)/2; %Half-size of the image
L = dx*C*2;
dx2 = dx*2*C/64;
aATF = abs(ATF); 
[x1,y1] = meshgrid(-L/2:dx:L/2-dx);
[x2,y2] = meshgrid(-L/2:dx2*8:L/2-dx);
typep = 'upsmp';
switch typep
    case 'upsmp'
        PSF = ifftshift(ifft2(ifftshift(aATF)));
        PSF2 = zeros(size(Img,1),size(Img,1));
        PSF2(C-32-1:C+31-1,C-32-1:C+31-1) = PSF*(C*2/64).^2;
        ATF2 = fftshift(fft2(fftshift(PSF2)));
        aATF2 = abs(ATF2);        
    case 'Hermite'
        %Hermite amplitude fitting
        aATF2 = zeros(size(Img));
        h = H64\aATF(:);
        aATF2_64 = aATF*0;
        aATF2_64(:) = H64*h;
        aATF2(:) = H256*h;
end




tmp_filter1 = aATF2 < max(aATF2(:))/exp(3);

%Zernike phase fitting
[Z1,~,is_in1] = zerfit( Img*0,x1,y1,0.14,Deg);
[Z2,~,is_in2] = zerfit( Img*0,x2,y2,0.14,Deg);

[ a2 ] = grad_svd( angle(ATF), 64, Z2, is_in2);        

% ATF(aATF2_64 < max(aATF2_64(:))/exp(3.5)) = 0;
a2(1:3) = 0;
Phase = zeros(size(Img));
Phase(is_in1) = Z1*a2;

%         Phase(aATF2 < max(aATF2(:))/exp(3.5)) = 0;

ATF2 = aATF2.*exp(1i*Phase);

%находим уровеь шума на изображении
fImg =  fftshift(fft2(fftshift(Img))) ;
noise = mean( abs(fImg(noise_area)) )/5;

afImg = imgaussfilt(abs(fImg),7);

SNR = (abs(afImg)-noise)./noise;

filtr = (1./(1 + 1./( SNR ) ) );
tmp_filter2 = filtr < max(filtr(:))/29;

% tmp_filter = tmp_filter2 | tmp_filter1; Убрал пока
tmp_filter = tmp_filter2; 
H2 = abs(ATF2).^2;

Dc = (1./ATF2).*(1./(1 + 1./( H2.*SNR.^2 ) ) ); %вот тут проблемы
Dc(tmp_filter) = 0;
Dc(isnan(Dc)) = 0;
Dc = Dc.*imgaussfilt(double(~tmp_filter1),25);
end

