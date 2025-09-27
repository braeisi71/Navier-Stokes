function y1 = f(X,t,type)
%KG Summary of this function goes here
%   Detailed explanation goes here
global v w
x = X(:,1);
y = X(:,2);
switch type
    case 'u1'
        y1 = -exp((-2*pi^2*t)*v).*sin(pi*x).*cos(pi*y);
    case 'u2'
        y1 = exp((-2*pi^2*t)*v).*cos(pi*x).*sin(pi*y);
    case 'u'
        u1 = t.^3.*y;
        u2 = t.^2.*x;
        y1 = [u1,u2];
end
end

