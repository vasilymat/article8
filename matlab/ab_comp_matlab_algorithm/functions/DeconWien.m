function [Dc] = DeconWien(ATF,Img,noise_area)
%Деконволюция с использованием фильтра виннера
%  

typep = 'interp';

switch typep
    case 'upsmp'
        PSF = ifftshift(ifft2(ifftshift(ATF)));
        PSF2 = zeros(256,256);
        PSF2(129-32-1:129+31-1,129-32-1:129+31-1) = PSF*(256/64).^2;
        ATF2 = fftshift(fft2(fftshift(PSF2)));
        aATF2 = abs(ATF2);
        % aATF2 = imgaussfilt(abs(ATF2),0);
        tmp_filter1 = aATF2 < max(aATF2(:))/exp(1)^4;
        ATF2 = aATF2.*exp(1j*angle(ATF2)); %Сгладил амплитуду, чтобы не было выбросов
    case 'interp'
        [X,Y] = meshgrid(1:4:256);
        [Xq,Yq] = meshgrid(1:1:256);
        ATF2 = interp2(X,Y,ATF,Xq,Yq,'makima');
        aATF2 = abs(ATF2);
        tmp_filter1 = aATF2 < max(aATF2(:))/exp(1)^6;
    case 'zerfit'
        dx = 300/256*10^-3;
        L = dx*256;
        [x1,y1] = meshgrid(-L/2:dx:L/2-dx);
        [x2,y2] = meshgrid(-L/2:dx*4:L/2-dx);
        [~,a2,~] = zerfit( ATF,x2,y2,0.14,12);
        [Z1,~,is_in1] = zerfit( Img*0,x1,y1,0.14,12);
        ATF2 = zeros(size(Img));
        ATF2(is_in1) = Z1*a2;
        aATF2 = abs(ATF2);
        tmp_filter1 = aATF2 < max(aATF2(:))/exp(1)^4;
        tmp_filter1 = tmp_filter1*0;
end
        
        
%находим уровеь шума на изображении
fImg =  fftshift(fft2(fftshift(Img))) ;
noise = mean( abs(fImg(noise_area)) );

afImg = imgaussfilt(abs(fImg),2);
SNR = (afImg-noise)./noise;

filtr = (1./(1 + 1./( SNR ) ) );
tmp_filter2 = filtr < max(filtr(:))/29;

tmp_filter = tmp_filter2 | tmp_filter1; 
H2 = abs(ATF2).^2;

Dc = (1./ATF2).*(1./(1 + 1./( H2.*SNR.^2 ) ) );
Dc(tmp_filter) = 0;
Dc(isnan(Dc)) = 0;
end

