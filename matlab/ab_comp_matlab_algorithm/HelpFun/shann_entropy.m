function [ Sh ] = shann_entropy( E )
%Получает массив и считает энтропию изображения

type_entropy = 1;

if size(E,3) > 1
    type_entropy = 5;
    [xdim,ydim,zdim] = size(E);
end


switch type_entropy
    case 1
        aT =(abs(E)).^2;
        aT = aT./sum(sum(aT))+10^-150; %чтобы не было деления на ноль. Тут м.б. ошибка
        Sh = -1*sum(sum(aT.*log2(aT)));
    
    case 2
        I = (abs( E./sum(sum(abs(E))).^2 ) );
        Sh = -1*sum(sum(I.*log2(I)));
    
    case 3
        E = abs(E);
        E = abs(E)/max(max(E))*255;
        %E = uint8(E);
        hE = hist(E,255);
        hE = sum(hE,2);
        [x,y] = size(E);
        hE = hE/x/y+10^-15;
        Sh = -sum(hE.*log2(hE));
    case 4
        aT =(abs(E)).^2;
        aT = aT./repmat( sum(sum(aT)), [xdim, ydim, 1] )+10^-150;
        Sh = -1*sum(sum(sum(aT.*log2(aT))));
    case 5
        aT =(sum(abs(E),3)).^2;
        aT = aT./sum(sum(aT))+10^-150; %чтобы не было деления на ноль. Тут м.б. ошибка
        Sh = -1*sum(sum(aT.*log2(aT)));        
end

end

