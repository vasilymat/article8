%������ �������� ��������
function  f = point_source(A,sig,k,x,y,xo,yo)
%������ �������
%F=(A/(sig*sqrt(2*pi)))*( exp( -1*( (x-xo).^2+(y-yo).^2 )/(2*sig^2) ) );%.*exp( 1i*k*sqrt( (x-xo).^2+(y-yo).^2 ) ) );

%����� �������
sig=2*max(max(x))/length(x)/1.5;
F=A*exp( -1*( (x-xo).^2+(y-yo).^2 )/(2*sig^2) ).*exp( -1i*( (x-xo).^2+(y-yo).^2 )/(2*sig^2)*0 );
%������� �� 
f=F;
end

