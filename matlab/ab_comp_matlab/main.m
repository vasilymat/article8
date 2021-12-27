% Program for finding and correction of complex image OTF usinf PGA method
% 

%Global constants and field prep
%%
addpath('C:\MatkivskiyV\science\article8\ab_comp_matlab\functions\')
addpath('C:\MatkivskiyV\science\article8\ab_comp_matlab\HelpFun\','-end')
%%
N = 256;
Global_constants = zeros(10,1);
lam = 0.85*10^-3; Global_constants(1)=lam; %Длина волны
k = 2*pi/(lam); Global_constants(2)=k; % Wave number in invers millimeters
dx = 300/256*10^-3; Global_constants(4)=dx; % it is pixel size in mm
L = dx*N; Global_constants(5)=L; % Length of image area

[x,y] = meshgrid(-L/2:dx:L/2-dx);
x = x + dx/2 ; % Это нужно, чтобы отцентрировать is_in_circle
y = y + dx/2;
kz = rphase3(dx,L,k);

Generate_HermPol = true; %Calculate Gauss-Hermite polynomials?
her_N = 6;    % The degree of Herm. pol.
her_R = 0.07; % The radius of Herm. pol.
fieldsprep;

radius = 0.14;
% This radius is determinated empirically. noise_area is area over which
% the noise is calculated
is_in_circle = circ(sqrt(x.^2+y.^2)/radius);
noise_area = circ(sqrt(x.^2+y.^2)/0.11) & ~circ(sqrt(x.^2+y.^2)/0.09);

Global_constants(10) = radius;
Global_constants(6) = 127; %xc
Global_constants(7) = 127; %yc

% tempz is used for sub-images that should not intersect with image borders
tempz = false(N,N);
tempz(33:223,33:223) = true;
clear Sh0;
%%

Q = [16,10,20,16];
N2 = size(Q,2);
Gst = repmat(E*0,[1 1 N2+1]);
Sh0(1) = shann_entropy(E(tempz));
u = fftshift(fft2(fftshift(E)));
Ev = E.*tempz;

fprintf('Initial image entropy = %g \n',Sh0(1));
for t2 = 1:1

    fprintf('circle %g \n',t2);
    fprintf('window = %g \n',Q(t2));

    [ATF] = PGA_find4_fft( Ev,Q(t2),tempz,radius,Global_constants,...
                                      8,radius,0, is_in_circle, 1, E);        
    sig = 0.15.^2;
    gauss = exp(-(x.^2+y.^2)/sig);

    Dc = DeconWienH(ATF,E,noise_area,H64,H256);
    u2 = u.*Dc.*gauss;
    u2 = u2*sqrt( sum(abs(u(:)).^2)/sum(abs(u2(:)).^2) );
    Ev = ifftshift( ifft2( ifftshift(u2)  ) );
%         u = u2;
    Sh0(t2+1) = shann_entropy(Ev(tempz));

    fprintf('Entropy inside the loop = %g \n',Sh0(t2+1));
end
[~,smin] = min(Sh0);
Gs = sum(Gst(:,:,1:smin),3);
Gk = ifft2( u.*ifftshift( exp(-1i*Gs)) );

figure(2); imagesc(abs(Ev)); colorbar; colormap(jet)
figure(3); imagesc(abs(u2)); colorbar; colormap(jet)
    



