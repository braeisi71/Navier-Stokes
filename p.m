function y1 = p(X,t)
%KG Summary of this function goes here
%   Detailed explanation goes here
global v w
x = X(:,1);
y = X(:,2);
y1 = .25*exp((-4*pi^2*t)*v).*(cos(2*pi*x)+cos(2*pi*y));
end

