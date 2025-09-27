function A = distance_matrix(X,Y)
[n,dim]=size(X);
m = size(Y,1);
A = zeros(n,m);
for d=1:dim
A = A + (repmat(X(:,d),1,m)-repmat(Y(:,d),1,n)').^2;
end
A = sqrt(A);