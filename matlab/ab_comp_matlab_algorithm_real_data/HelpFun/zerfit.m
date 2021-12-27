function [Z,a,is_in_circle] = zerfit( Phase,x,y,R,Pol_deg )
%zerfit( Массив для апроксимации,x,y,радиус апертуры,степень полиномов (до 7ми включительно)
[theta,r]=cart2pol(x,y);
N=[]; M=[];
for n=0:Pol_deg
    N = [N n*ones(1,n+1)];
    M = [M -n:2:n];
end
r=r/R;
is_in_circle = ( r <= 1 );
Z = zernfun(N,M,r(is_in_circle),theta(is_in_circle));
a = Z\Phase(is_in_circle);
end

