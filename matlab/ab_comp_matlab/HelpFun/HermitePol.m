function [HermArr] = HermitePol(x,y,N,w)
%Generating an 2D array of Hermite polynomials
%size - 1D shape of array
%N - maximum degree of a monomial.
HermArr =  zeros(length(x(:)),N);
tmpArrX = zeros(size(x,1),size(x,1),N);
tmpArrY = zeros(size(x,1),size(x,1),N);

s2w = sqrt(2)/w;
w2 = 1/w^2;
expprod = exp(-w2*x.^2).*exp(-w2*y.^2);


% parfor n = 1:N
%     tmp = hermiteH(n-1,s2w*x);
%     HermArr(:,n) = tmp(:);
% end


parfor n = 1:N   
    tmp = hermiteH(n-1,s2w*x);
    tmpArrX(:,:,n) = tmp;
    tmpArrY(:,:,n) = tmp';
end


for n1 = 1:N
    for n2 = 1:N
        tmp = tmpArrX(:,:,n1).*tmpArrY(:,:,n2).*expprod;
        HermArr(:,n2 + (n1-1)*N) = tmp(:);
    end
end


% for n1 = 1:N
%     for n2 = 1:N
%         tmp = hermiteH(n1-1,s2w*x).*hermiteH(n2-1,s2w*y).*expprod;
%         HermArr(:,n2 + (n1-1)*N) = tmp(:);
%     end
% end

end

