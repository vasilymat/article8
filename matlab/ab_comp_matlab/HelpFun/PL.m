function f=PL(dz,kz,T)
%������� ������������ ������������ ����

%� �������������� ������� ��������
f=ifftshift( ifft2( fft2(fftshift(T)).*fftshift(exp(1i*kz*dz)) ) );
%�����

% %� �������������� ������� ������ �� �������
% dx=5.5*10^-3;
% lam=0.0008477;
% k=2*pi/lam; 
% h=1/(1i*lam*dz)*exp(1i*k*kz/(2*dz));
% H=fft2(fftshift(h))*dx^2;
% U1=fft2(fftshift(T));
% U2=U1.*H;
% f=ifftshift(ifft2(U2));
% %�����
end
