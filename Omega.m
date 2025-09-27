function Points = Omega(a,b,c,d,h,dist_type,domain_type)
switch nargin 
    case 6
        domain_type = 'regular';
    case 5
        domain_type = 1';
        dist_type = 'Uniform';
end
switch domain_type
    case 1 % Rectangular domain
        yh=a+h:h:b-h;
        yv=c+h:h:d-h;
        n1=size(yh,2);
        n2=size(yv,2);
        GL=[a*ones(n2,1),yv'];
        GR=[b*ones(n2,1),yv'];
        GU=[yh',d*ones(n1,1)];
        GD=[yh',c*ones(n1,1)];
        switch dist_type
            case 'Halton'
                K  = floor(((b-a)/h-1)*((d-c)/h-1));
                HaltonPoints = Halton(K,a+h/2,b-h/2,c+h/2,d-h/2);
                Xi = HaltonPoints(:,1);
                Yi = HaltonPoints(:,2);
                Y = [Xi,Yi];
            case 'Uniform'
                [Xi,Yi]= meshgrid(yh',yv');
                Y=[Xi(:), Yi(:)];
        end
        Points{1,1} = [Y(:,1),Y(:,2)];
        Points{1,2} = 'Interior points';
        Points{2,1}=[GL;GR;GU;GD;a,c;b,c;a,d;b,d];
        Points{2,2} = 'Boundary points';
        Points{3,1} = GL;
        Points{3,2} = 'Left boundary';
        Points{4,1} = GR;
        Points{4,2} = 'Right boundary';
        Points{5,1} = GU;
        Points{5,2} = 'Upper boundary';
        Points{6,1} = GD;
        Points{6,2} = 'Down boundary';
        Points{7,1} = [a,c;b,c;a,d;b,d];
        Points{7,2} = 'Vertices';
end
function Points = Halton(N,a,b,c,e,d)
switch nargin
    case 0
        N = 1000;
        d = 2;
        a = 0; b = 1; c=0; e=0;
    case 1
        d = 2;
        a = 0; b = 1; c=0; e=1;
    case 5
        d=2;
end
p = haltonset(d,'Skip',1e3,'Leap',1e2);
X = net(p,N);
x=(b-a)*X(:,1)+a;
y=(c-e)*X(:,2)+e;
Points = [x,y];



