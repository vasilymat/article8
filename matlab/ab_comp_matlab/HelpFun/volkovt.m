function [ fi ] =volkovt( gamma2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

linhas=size(gamma2,1);
colunas=size(gamma2,2);

image=gamma2;
image=padarray(image,[abs(linhas-colunas),abs(linhas-colunas)],'replicate','post');%Faz padarray se a imagem nao for quadrada
if mod(size(image,1),2)~=0
image=padarray(image,[1,0],'replicate','post');%Faz padarray se as dimensoes nao forem par
end
if mod(size(image,2),2)~=0
image=padarray(image,[0,1],'replicate','post');%Faz padarray se as dimensoes nao forem par
end


[sy,sx]=size(image);
p=linspace(-sx/2,sx/2-1,sx);
q=linspace(-sy/2,sy/2-1,sy);
[p,q]=meshgrid(p,q);
N=sx;

% figure;
% imshow(image,[]);



Gpsix=cos(image).*ifft2(fftshift((p*2*pi).*fftshift(fft2(sin(image)))))...
    -sin(image).*ifft2(fftshift((p*2*pi).*fftshift(fft2(cos(image)))));
Gpsiy=cos(image).*ifft2(fftshift((q*2*pi).*fftshift(fft2(sin(image)))))...
    -sin(image).*ifft2(fftshift((q*2*pi).*fftshift(fft2(cos(image)))));

Gfix=ifft2(fftshift((p*2*pi).*fftshift(fft2(image))));
Gfiy=ifft2(fftshift((q*2*pi).*fftshift(fft2(image))));

Gkx=(Gpsix-Gfix)/2/pi;
Gky=(Gpsiy-Gfiy)/2/pi;

p(sy/2+1,sx/2+1)=1;
b=(fftshift(fft2(Gkx)).*p+fftshift(fft2(Gky)).*q)./(p.^2+q.^2);

k=(ifft2(fftshift(b)));

fi=double(real(k)+image);
% fi=gather(fi);

% 
% figure;
% imshow(fi,[]);
% 
% figure;
% surfc(fi); shading interp;
fi=fi(1:size(gamma2,1),1:size(gamma2,2));

end