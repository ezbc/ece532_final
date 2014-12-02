function [ xhat ] = Lasso2( A,lambda,b,maxIter,eps )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

xhat = zeros(size(A,2),1);
size(lambda*(A'*A) )
alpha = 1./(lambda*(A'*A)); %Something seems wrong here.
delta = 10;
iter = 0;

while( iter < maxIter && delta >eps)
    y = xhat + alpha*A'*(b-A*xhat);
    xnext = sign(y) * max(abs(y) - alpha*lambda,zeros(size(A,2),1));
    iter = iter+ 1;
    delta = norm(xnext - xhat);
    xhat = xnext;
end
end

