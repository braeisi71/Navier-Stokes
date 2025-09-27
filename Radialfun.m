function A = Radialfun(Xe,X,rbftype,delta,op)
%توابع پایه شعاعی را حساب میکنه
% rbfpar for wendland functions is Kind of wendland function
% rbfscale is the scale of trial function
% operator is the order of partial derivative of function
%rbf_kind = 'Wendland','Gaussian', 'MQ', 'IMQ', 'Polyharmonic'
global rbfpar rbfscale
n = size(Xe,1);
[m,dim] = size(X);
r  = distance_matrix(Xe,X)+eps;
x = diff_matrix(Xe,X,'x');
y = diff_matrix(Xe,X,'y');
%%
switch rbftype
    case 'Wendland32'% phi_3,2   in C^4
        e = 1/delta;
        y1 = max(1-e*r,0);
        switch op
            case '0'
                A = y1.^6.*(35*(e*r).^2+18*(e*r)+3);
            case 'x'
                A = -56*e^2*y1.^5.*(5*e*r+1).*x;
            case 'y'
                A = -56*e^2*y1.^5.*(5*e*r+1).*y;
            case 'L'
                phi_rDivr = -56*e^2.*y1.^5.*(5*e*r+1);
                phi_rr = 56*e^2*y1.^4.*(35*(e*r).^2-4*e*r-1);
                A = phi_rr + phi_rDivr;
        end
    case 'MQ'
        e = delta;
        switch op
            case '0'
                A = sqrt(1+(e*r).^2);
            case 'x'
                phi_rDivr = (e^2)./sqrt(1+(e*r).^2);
                A = x.*phi_rDivr;
            case 'y'
                phi_rDivr = (e^2)./sqrt(1+(e*r).^2);
                A = y.*phi_rDivr;
            case 'L'
                phi_rDivr = (e^2)./sqrt(1+(e*r).^2);
                phi_rr = e^2./(1+(e*r).^2).^(3/2);
                A = phi_rDivr + phi_rr;
        end 
end


