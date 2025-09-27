function [Xall,Xin,Xbd] = points_sq_unity(a,b,c,d,h,type)
%   Ommega: the points in the region
%   Gamma1: the points on the inner boundary of the region
%   Gamma2= the points on the outer boundary
%   N= number of total points
switch type
    case 'R'
        yh=a-eps:h:b+eps;
        yv=c-eps:h:d+eps;
        [Y1,Y2]= meshgrid(yh',yv');
        Y=[Y1(:), Y2(:)];
        Xin=[Y(:,1),Y(:,2)];
        yh=a:h:b;
        yv=c:h:d;
        [Y1,Y2]= meshgrid(yh',yv');
        Xall=[Y1(:), Y2(:)];
        n1=size(yh,2);
        n2=size(yv,2);
        GL=[a*ones(n2,1),yv'];
        GR=[b*ones(n2,1),yv'];
        GU=[yh(2:end-1)',d*ones(n1-2,1)];
        GD=[yh(2:end-1)',c*ones(n1-2,1)];
        Xbd=[GL;GR;GU;GD];
end
