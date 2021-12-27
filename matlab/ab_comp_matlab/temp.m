if 1>0.5
    .5
end
if 1>.3
    .3
end


% %Попытаемся понять, как разложить имеющуюся ATF по полиномам эрмита
% her_N = 6;
% her_R = 0.07;
% dx2t = dx*N/64;
% if Generate_HermPol == true
%     [x2,y2] = meshgrid(-L/2:dx2t:L/2-dx2t);
%     x2 = x2 + dx/2; % Это нужно, чтобы отцентрировать is_in_circle
%     y2 = y2 + dx/2;    
%     H64 = HermitePol(x2,y2,her_N,her_R);
% %     H256 = HermitePol(x,y,her_N,her_R);
% end
% 
% aATF2 = zeros(size(x2));
% h = H64\aATF(:);
% aATF2(:) = H64*h;
% figure(1); imagesc(aATF2); colorbar; colormap(jet)
% %%
% aATF3 = zeros(size(x));
% aATF3(:) = H256*h;
% 
% figure(2); imagesc(aATF3); colorbar; colormap(jet)