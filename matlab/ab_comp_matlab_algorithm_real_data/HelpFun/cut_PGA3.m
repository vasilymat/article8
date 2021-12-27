function [ M,P2,cx,cy ] = cut_PGA3( P2,sT2,exclude_new,da )
%Вырезаем в общем

da = da*2;
M = max(sT2(P2));
[cx,cy] = find( sT2 == M);
P2(cx-da:cx+da-1,cy-da:cy+da-1) =...
    P2(cx-da:cx+da-1,cy-da:cy+da-1) & ~exclude_new;


end

