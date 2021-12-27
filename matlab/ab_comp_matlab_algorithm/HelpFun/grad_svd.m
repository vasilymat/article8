function [ a3 ] = grad_svd( aP, l2, Z, is_in_circle)
%Функция ищет коэффициенты по градиентам фазы Se
%в разложении использует и константу с наклоном
Se = zeros(l2,l2);
size_z = size(Z,2);
a = zeros(size_z,1);

aP(~is_in_circle) = nan;

aP_x = diff(aP,1,1);
aP_y = diff(aP,1,2);

is_in_circle_x = ~ isnan(aP_x); %находим, где не nan
is_in_circle_y = ~ isnan(aP_y);

if_overpi = aP_x >= pi;
aP_x(if_overpi) = aP_x(if_overpi) - 2*pi;
if_lesspi = aP_x <= -pi;
aP_x(if_lesspi) = aP_x(if_lesspi) + 2*pi;

if_overpi = aP_y >= pi;
aP_y(if_overpi) = aP_y(if_overpi) - 2*pi;
if_lesspi = aP_y <= -pi;
aP_y(if_lesspi) = aP_y(if_lesspi) + 2*pi;

aP_x(~is_in_circle_x) = 0;
aP_y(~is_in_circle_y) = 0;


       
% тут мы определяем количество и место ненулевых элементов после
% дифференцирования
gr_numx = length(aP_x(is_in_circle_x));
gr_numy = length(aP_x(is_in_circle_y));
M = zeros(gr_numx+gr_numy,size_z-1);

s = zeros(gr_numx+gr_numy,1);
s(1:gr_numx) = aP_x(is_in_circle_x);
s(gr_numx+1:end) = aP_y(is_in_circle_y);

    for ii = 2:size_z
        a(ii) = 1;
        Se(is_in_circle) = Z*a;
  
        Se_x = diff(Se,1,1);
        Se_y = diff(Se,1,2);        
      
        M(1:gr_numx,ii-1) = Se_x(is_in_circle_x);
        M(gr_numx+1:end,ii-1) = Se_y(is_in_circle_y);           

        a(ii) = 0;
    end
    
[U, D, V] = svd(M,0);
a3(2:size_z) = (V)*pinv(D)*(U')*s;

a3 = a3';


end

