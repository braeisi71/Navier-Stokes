function y1 = g(X,t,type)
%KG Summary of this function goes here
%   Detailed explanation goes here
global v w
x = X(:,1);
y = X(:,2);
switch type
    case 'u1'
        y1 = (exp(-2*pi^2*t*v)).^2 .* sin(pi*x) .* cos(pi*y).^2 .* pi .* cos(pi*x) + ...
            (exp(-2*pi^2*t*v)).^2 .* cos(pi*x) .* sin(pi*y).^2 .* sin(pi*x) .* pi - ...
            (1/2) * exp(-4*pi^2*t*v) .* pi .* sin(2*pi*x);
    case 'u2'
        y1 = (exp(-2*pi^2*t*v)).^2 .* sin(pi*x).^2 .* cos(pi*y) .* pi .* sin(pi*y) + ...
            (exp(-2*pi^2*t*v)).^2 .* cos(pi*x).^2 .* sin(pi*y) .* pi .* cos(pi*y) - ...
            (1/2) * exp(-4*pi^2*t*v) .* pi .* sin(2*pi*y);
end

