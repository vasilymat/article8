%Линза с жестким диаметром 40
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Диаметр40\Aber21MZ31_40b.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Диаметр40\Aber21MZ00_40.mat')
% load('C:\MatkivskiyV\!Актуальная наука\статья 8\matlab\Диаметр40\Aber21MZ20_40.mat')


%==real aberrations from setup==
%read no-aberrated data
D1 = fopen('C:\MatkivskiyV\science\article8\ab_comp_matlab\data\data_from_setup\USAF_4.dat');
A = fread(D1,1024*256*256*2,'single');
fclose(D1);
B = A(1:2:end) + 1j*A(2:2:end);
C = reshape(B,[256,256,1024]);
E0 = squeeze(sum(C(:,105,1:2:end),2));
E = zeros(size(E0),'single');
E(73:203,257:472) = E0(73:203,257:472);
clear A B C E0
% centring
fE = FTS(E);
fE  = circshift(fE,[135 -20]);
fE2 = zeros(512,512);
fE2(129:384,:) = fE(:,:);
E =  circshift((FTS(fE2)),[43 111]);
E = E*256./max(abs(E(:)));

%read aberrated data
D1 = fopen('C:\MatkivskiyV\science\article8\ab_comp_matlab\data\data_from_setup\USAF_4_ab_2.dat');
D = zeros(1024,256,256,'single');
A = fread(D1,1024*256*256*2,'single');
fclose(D1);

B = A(1:2:end) + 1j*A(2:2:end);
C = reshape(B,[256,256,1024]);
E0a = squeeze(sum(C(:,108:112,1:2:end),2));
Ea = zeros(size(E0a),'single');
Ea(73:203,228:444) = E0a(73:203,228:444);
clear A B C E0a
fEa = FTS(Ea);
fEa  = circshift(fEa,[135 -20]);
fE2a = zeros(512,512);
fE2a(129:384,:) = fEa(:,:);
Ea2 =  circshift((FTS(fE2a)),[22 75]);
Ea = Ea2*256./max(abs(Ea2(:)));
clear Ea2

