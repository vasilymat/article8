function [ a3 ] = grad_svd2( aP_x, aP_y, l2, Z, is_in_circle, is_in_circle2)
%Функция ищет коэффициенты по градиентам фазы Se
%в разложении использует и константу с наклоном
%В отличии от предыдущей - получает полные градиенты сразу, а не вычисляет
%их внутри
Se = zeros(l2,l2);
size_z = size(Z,2);
a = zeros(size_z,1);

temp = zeros(64,64);
temp(~is_in_circle) = nan;

temp_x = diff(temp,1,1);
temp_y = diff(temp,1,2);

is_in_circle_x = ~ isnan(temp_x); %находим, где не nan
is_in_circle_y = ~ isnan(temp_y);


 
% figure(1); imagesc(aP_x.*is_in_circle_x); colorbar ; colormap(gray);
% figure(2); imagesc(aP_y.*is_in_circle_y); colorbar ; colormap(gray);
        
% тут мы определяем количество и место ненулевых элементов после
% дифференцирования
gr_num = length(aP_x(is_in_circle_x));
M = zeros(2*gr_num,size_z-1);

s = zeros(2*gr_num,1);
s(1:gr_num) = aP_x(is_in_circle_x);
s(gr_num+1:end) = aP_y(is_in_circle_y);

    for ii = 2:size_z
        a(ii) = 1;
        Se(is_in_circle) = Z*a;        
       
         Se(~is_in_circle) = nan;    
         Se_x = diff(Se,1,1);
         Se_y = diff(Se,1,2);
      
        M(1:gr_num,ii-1) = Se_x(is_in_circle_x);
        M(gr_num+1:end,ii-1) = Se_y(is_in_circle_y);           

        a(ii) = 0;
    end
    
[U, D, V] = svd(M,0);
a3(2:size_z) = (V)*pinv(D)*(U')*s;

%a3(2:size_z) = (M'*M)\M'*s;

a3 = a3';


end

