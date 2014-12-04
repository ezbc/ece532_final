function [ xhat ] = Lasso( A,lambda,b,maxIter,eps )
%Lasso computes the xhat that solves min norm2(A-X*b)- lamnda * norm(x).
% This done so that the resulting xhat is mostly zero and thus chooses the
% appropiate features that we need from a sparse amtrix A.
%
% More info at http://comisef.eu/files/wps042.pdf
%              http://statweb.stanford.edu/~tibs/lasso/lasso.pdf
%              https://www.ics.uci.edu/~welling/people/lboyles_thesis.pdf
%              (pages 10-11)

Xcurrent = ones(size(A,2) ,1);
XPrev = inf(size(A,2) ,1);
currentIter = 0;

while( currentIter < maxIter && norm(Xcurrent - XPrev) > eps )
    currentIter = currentIter + 1
    XPrev = Xcurrent;
    
    for j = 1:size(Xcurrent,1)
        
        bhat = sum(A*Xcurrent) -A(:,j)*Xcurrent(j) ;
        
        
        alpha(j) = 1/size(A,1) * sum(A(:,j)'*(b-bhat));
        if(alpha(j) > 0 && lambda <abs(alpha(j)) )
            S(j) = alpha(j) - lambda;
        elseif(alpha(j) < 0 && lambda <abs(alpha(j)) )
            S(j) = alpha(j) + lambda;
        else
            S(j) = 0;
        end
    end

    Xcurrent = S(:,:)./sum(A.^2);
    
    Xcurrent = Xcurrent';
    
    
    %removed the NaN value to zero as we divided by zero somewhere
    Xcurrent(find(isnan(Xcurrent))) = 0;
    
    
    %
    %     W = diag(abs(Xcurrent));
    %     [T,S,U] = svd(W)
    %     size(A'*A)
    %     size(W)
    %     size(A' * b)
    %
    %     Xcurrent = (A'*A + lambda * pinv(W))^-1 * A' * b;
end

xhat = Xcurrent;
end

