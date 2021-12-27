function [norm] = normf(E,E0)
% E is normalized function
% E0 is normilizing function (i.e. function, that norm is standart)
    norm = sqrt( sum(abs(E0(:)).^2)./sum(abs(E(:)).^2) );
end 
