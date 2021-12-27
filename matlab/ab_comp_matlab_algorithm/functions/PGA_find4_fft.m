function [ f ] = PGA_find4_fft( T,sq,tempz,Rex,Global_constants,...
    TypeWindow,radius, Z, is_in_circle, Deg, E2 )

% В этом файле делаем процедуру разбиение на локальные радиусы для eig
% radius - это текущий радиус
% R - минимальное расстояние между центрами квадратов
% is_in_circle - передает текущую полную апертуру
% prev_phase - фаза, найденная на предыдущем этапе
% Rex - радиус окружности для исключения

k=Global_constants(2);
dx=Global_constants(4); 
L=Global_constants(5);
Rf = Global_constants(10); %полный радиус
nsq = 15;%оценочно - максимальное количество точек
l = size(T,2);
f = zeros(size(T));
%the end of block of global constants

%формулы для собств. вектора
l2 = 64;
dx2 = dx*l/l2;
[x2,y2] = meshgrid(-L/2:dx2:L/2-dx2);

x2 = x2 + dx2/2; % Это нужно, чтобы отцентрировать is_in_circle
y2 = y2 + dx2/2;

tempz2 = tempz;
% tempz2 = false(1024,1024);
% tempz2(301:807,241:770) = true;

Ps2 = winmultp(0.1); %Задаем множитель, который будт определять ширину окна
    function Psm = winmultp(sig)        
        is_in_circle_t = circ(sqrt((x2-dx2/2).^2+(y2-dx2/2).^2)...
                                                            /(L/2*sq/l2));
        Psm = imgaussfilt(double(is_in_circle_t),sq*sig);
    end
 
%local constants

is_eig = true;
Npie = 1;   %Это на такое количество отрезков делится размерность
rg = 32/1;  %r_grid - задаем радиус субобласти для eig
rg_o = 1; %было 8 в обоих случаях   %rg_oder, tau_o - это для третьего случая (type_eigpiece = 3) 
tau_o = 1;  %- на сколько элементов делится апертура по tau и r.

type_eigpiece = 3;

% A=zeros(150,2);%Массив для настоящих максимумов, откуда будет вырезан квадратик
%the end of block of local constants
%------------------------------

da=fix(sq/2); %dELTa

% в этом куске формируются переменные для eig по частям

is_in_circle2 = circ(sqrt((x2).^2+(y2).^2)/Rf);
is_in_circle4 = circ(sqrt((x2).^2+(y2).^2)/(radius-0*dx2)*1);

[C, en, masl, Ntemp, Z2] = eigpiece(type_eigpiece);
    function [C, en, masl, Ntemp, Z2] = eigpiece(typep)
        exclude_is_in = circ(sqrt((x2).^2+(y2).^2)/Rex);
        is_in_circle3 = ~exclude_is_in & is_in_circle2; %это незадействованные ранее точки
%         C = zeros( length(x2(is_in_circle3)).^2,1); %это строка, cодержащая всевозможные строки изображения
%         en = zeros(Npie*Npie+1,1,'uint32');
        [Z2,~,~] = zerfit(x2,x2,y2,Rf,Deg); %тут мы сгенирировали фазу, которую
        
        
        switch typep 
          case 1
            masl = false(l2,l2,Npie*Npie);
            masch = zeros(l2/Npie,Npie,'int8');
            for ii2 = 1:Npie
                masch(:,ii2) = 64*(ii2-1)/Npie+1:64*ii2/Npie;
            end        
    
            for ii2 = 1:Npie    
                for jj = 1:Npie
                    Ntemp = (ii2-1)*Npie+jj;
                    masl(masch(:,ii2),masch(:,jj),Ntemp) = masl(...
                        masch(:,ii2),masch(:,jj),Ntemp) |...
                        is_in_circle3(masch(:,ii2),masch(:,jj));
        
                    %end - конец соответствующего участка
                    en(Ntemp+1) = size(x2(masl(:,:,Ntemp)),1);
         
                end
            end
            en = en.^2;
            en = cumsum(en);
            %конец eig по частям
          case 2 %!альтернативный eig по частям

            dots = is_in_circle2*0;
            sum_dots = dots;


            rs = fix(rg*sqrt(2)*0.99); % расстояние между центрами субобластей в sqrt2 больше радиуса
                           % сделано, чтобы не было пустых областей по
                           % диагонали

            deg = fix(l2/2/rg-10^-10); % это сколько вправо и влево точек сетки 0
                           % -10^-10  -для того, чтобы "вычесть границу"
                           % 32/32 должно быть меньше 1

            masl = false(l2,l2); %это битовый трехмерный массив с масками каждой из субобластей
            Ntemp = 0;
            for ii2 = -deg:deg
                for jj = -deg:deg
                    if (sqrt( (ii2*rs).^2 + (jj*rs).^2) <= l2/4) 
                        Ntemp = Ntemp + 1;
                        dots(33-ii2*rs,33-jj*rs) = true; %записываем 1 в центры областей
                        temp_circ = circ(sqrt((x2/dx2-ii2*rs).^2 ...
                             +(y2/dx2-jj*rs).^2)/rg);
                        temp_circ = temp_circ & is_in_circle2; 
                        masl(:,:,Ntemp) = temp_circ;
                        length(x2(temp_circ));
                        sum_dots = sum_dots + temp_circ; % это массив, на которую 
                                                         % потом отнормируем сум eig
                        en(Ntemp+1) = size(x2(temp_circ),1);
                    end
                end
            end
            sum_dots = sum_dots - 0.5; % тут избавляемся от деления на ноль
            sum_dots = abs(sum_dots);
            sum_dots = sum_dots + 0.5;

            en = en.^2;
            en = cumsum(en);
            %!конец альтернативный eig по частям
            
          case 3
            masl = false(l2,l2); %это битовый трехмерный массив с масками
                                 %каждой из субобластей
            [tau,rn] = cart2pol(x2,y2);
            
            rt = radius/rg_o;     %радиус центральной части
            dtau = 2*pi/tau_o; %шаг по tau
            
            temp_circ = rn <= rt;
            masl(:,:,1) = temp_circ;
            Ntemp = 1;
            en(Ntemp+1) = size(x2(temp_circ),1);
            
            for ii2 = 1:rg_o-1
                temp_circ2 = (rn <= radius*(1/rg_o+ii2/rg_o)) & ~temp_circ;
                for jj = 1:tau_o
                   Ntemp = Ntemp + 1;
                   temp_circ_tau = (tau >= -pi+(jj-1)*dtau) &...
                       (tau < -pi+jj*dtau);
                   masl(:,:,Ntemp) = temp_circ_tau & temp_circ2;
                   %figure(1); imagesc(masl(:,:,Ntemp)); colorbar ; colormap(gray);
                   
                   en(Ntemp+1) = size(x2(masl(:,:,Ntemp)),1);
                   
                end
                temp_circ = temp_circ2 | temp_circ;
            end
            
            sum_dots = temp_circ*0+1;
            en = en.^2;
            en = cumsum(en);
            C = zeros( en(end),1);
        case 4
            masl = false(l2,l2);
            Ntemp = 0;
            temp_circ = false(l2,l2);
            [tau,rn] = cart2pol(x2,y2);
%             rt = radius/rg_o; %Главный субрадиус
            for ii2 = 0:rg_o-1
                temp_circ2 = (rn <= radius*(1+ii2)/rg_o) & ~temp_circ;
                tau_o = 1 + ii2*2;
                dtau = 2*pi/tau_o; %шаг по tau
                for jj = 1:tau_o
                    Ntemp = Ntemp + 1;
                    temp_circ_tau = (tau >= -pi+(jj-1)*dtau) &...
                       (tau < -pi+jj*dtau);
                    masl(:,:,Ntemp) = temp_circ_tau & temp_circ2;
                    en(Ntemp+1) = size(x2(masl(:,:,Ntemp)),1);
                end
                temp_circ = temp_circ2 | temp_circ;
            end
            en = en.^2;
            en = cumsum(en);
            C = zeros( en(end),1);
            
        end
        
        
    end

dat = da;
da = fix(da/1);
[x3,y3] = meshgrid(-1:1/(2*da):1-1/(2*da));
x3 = x3 + 1/(4*da); y3 = y3 + 1/(4*da);
exclude_new = circ(x3.^2+y3.^2);
da = dat;

% A = A*0;   
P2 = tempz2;
if sq > 50
    sT = abs(E2);
else
    sT = abs(T);
end
stop_while = true;
G = zeros(l2,l2,250);
G1 = G;
sh = 0;
da2 = l2/2;

N = 1;
while stop_while && (10 > N) && any(P2(:)) 

    [mas_max(N),P2,cx,cy] = cut_PGA3(P2,sT,exclude_new,fix(da/1));

    if mas_max(1)/mas_max(N) > 1.5
        stop_while = false;
    end
    A(N,1) = cx(1);
    A(N,2) = cy(1);
    %P1 - массив комплексных чисел, представляет собой вырезанный квадратик
    %P2 - также массив комплексных чисел, содержит всевозможные вырезы

    %Квадратный вырез для окна l2
    %Умножаем на простой Гаусс по ширине окна
%     [cx,cy] %Тут что-то нужно придумывать!(((
    
    cxm = cx+sh-da2   ;
    cxp = cx+sh+da2-1 ;
    cym = cy+sh-da2   ;
    cyp = cy+sh+da2-1 ;
    sx1 = 1;
    sx2 = 64;
    sy1 = 1;
    sy2 = 64;
    
    if cx+da2-1 > 256
        cxp = 256;
        sx2 = 64 - (cx+da2-1-256);
    end
    if cx-da2 < 1
        cxm = 1;
        sx1 = da2-cx+2;
    end
    if cy+da2-1 > 256
        cyp = 256;
        sy2 = 64 - (cy+da2-1-256);
    end
    if cy-da2 < 1
        cym = 1;
        sy1 = da2-cy+2;
    end
           
    S3 = zeros(64,64);    
    S3(sx1:sx2,sy1:sy2) = T(cxm:cxp,cym:cyp);        
    G(:,:,N) = fftshift(fft2(fftshift(S3.*Ps2)));
    N = N + 1;
end


fprintf('Количество точек = %g \n',N);

sum_um = zeros(size(G(:,:,1)));
sum_vm = zeros(size(G(:,:,1)));

if is_eig == true    

    for kk = 1:N-1
%         Gt = exp(1i*angle(G(:,:,kk)));
        Gt = G(:,:,kk);
        for ii = 1:Ntemp        
            xv = Gt(masl(:,:,ii));           
            Ct = xv*xv';
            C(en(ii)+1:en(ii+1)) = C(en(ii)+1:en(ii+1)) + Ct(:);
        end

    end
  
else
    for kk = 1:N-1
        Gt = exp(1i*angle(G(:,:,kk)));
        sum_um = sum_um + Gt.*conj( circshift(Gt,[1 0]) );
        sum_vm = sum_vm + Gt.*conj( circshift(Gt,[0 1]) );
    end
    aP_x = angle(sum_um).*is_in_circle2;
    aP_y = angle(sum_vm).*is_in_circle2;

    f3 = r_hand3(aP_x(1:end-1,:),aP_y(:,1:end-1));
end


P = zeros(l2,l2);
    
en2 = diff(en);
en2 = int32(sqrt(single(en2)));


if is_eig == true 
    for ii = 1:Ntemp %!!! ВОт тут происходит обратный сбор
        dsize = en2(ii);
        Ct = zeros(dsize,dsize);
        Ct(true(dsize,dsize)) = C(en(ii)+1:en(ii+1));
        [V,~] = eig(Ct);
        P(masl(:,:,ii)) = V(:,end);
        f = P;
    end
    
else
    a_end = f3(is_in_circle2)\Z2;
    a_end(1:3) = 0;
    f(is_in_circle) = Z*a_end';
    
end


end




