function [Dc] = WinF(ATF2,Img,noise_area)
% Creating deconvolution multiplier from Amplitude transfer fanction 
% ATА2 - full grid ATF
% Img - aberrated image
% noise_area - disk bool area, that defines noise in image spectr
aATF2 = abs(ATF2);
tmp_filter1 = aATF2 < max(aATF2(:))/exp(3);

fImg =  fftshift(fft2(fftshift(Img))) ;
noise = mean( abs(fImg(noise_area)) )/15;

afImg = imgaussfilt(abs(fImg),7);
% ht = H256\afImg(:);
% afImg(:) = H256*ht;

SNR = (abs(afImg)-noise)./noise;

filtr = (1./(1 + 1./( SNR ) ) );
tmp_filter2 = filtr < max(filtr(:))/29;

% tmp_filter = tmp_filter2 | tmp_filter1; Убрал пока
tmp_filter = tmp_filter2; 
H2 = aATF2.^2;

Dc = (1./ATF2).*(1./(1 + 1./( H2.*SNR.^2 ) ) ); %вот тут проблемы
Dc(tmp_filter) = 0;
Dc(isnan(Dc)) = 0;
Dc = Dc.*imgaussfilt(double(~tmp_filter1),25);
end

