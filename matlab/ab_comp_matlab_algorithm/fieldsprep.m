% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Aber21MZ31.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Aber21MZ0.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Aber21MZ20_120.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Aber21MZ20_40.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Aber21MZ31_100.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Aber21MZ20_900.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Aber21MZ31_900.mat')

%Линза с жестким диаметром 40
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Диаметр40\Aber21MZ31_40b.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Диаметр40\Aber21MZ00_40.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Диаметр40\Aber21MZ20_40.mat')


%==real aberrations==
%10 doths without abb
load('C:\MatkivskiyV\science\article8\ab_comp_matlab\data\modelling_real_abb\OCT10_Bez.mat')
Eideal = L2(:,:,51);
%10 doths without low abb
load('C:\MatkivskiyV\science\article8\ab_comp_matlab\data\modelling_real_abb\OCT10_Se_lNa256.mat')
E = L2(:,:,51);

%Линза с жестким диаметром 40 и без аберрация на проход туда
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\файлы с хорошим пучком\Aber21MZ00_40_bez.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\файлы с хорошим пучком\Aber21MZ31_40b_bez.mat')

%Большая линза
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\файлы через большую линзу\Aber21GZ00_40.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\файлы через большую линзу\Aber21GZ31_40.mat')


% type_prep = 'Aber21MZ31_900.mat';

% T = L2(:,:,114);
% 
% switch type_prep
%     case 'Aber21MZ31_900.mat'
%         fE = fftshift(fft2(T));
%         fE = circshift(fE, [-60,0]);
%         E = ifft2(ifftshift(fE));
%     case 'else'
%         fT = fft2(T);
%         fT = circshift(fT, [-60,0]);
%         E = ifft2((fT));
% end

% %noise insertion
% tmp = randn([256,256])*30;
% E0 = E;
% E = E + tmp;

% figure(1); imagesc(abs(E) ); colorbar; colormap(jet);


% dx2t = dx*N/64;
% if Generate_HermPol == true
%     [x2,y2] = meshgrid(-L/2:dx2t:L/2-dx2t);
%     x2 = x2 + dx/2; % Это нужно, чтобы отцентрировать is_in_circle
%     y2 = y2 + dx/2;
%     H256 = HermitePol(x,y,her_N,her_R);
%     H64 = HermitePol(x2,y2,her_N,her_R);
% end

