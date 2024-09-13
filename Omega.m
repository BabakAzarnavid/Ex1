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
    case 2 % Eliptical domain
        switch dist_type
            case 'Halton'
                N  = floor(((b-a)/h)*((d-c)/h));
                HaltonPoints = Halton(N,-1,1,-1,1);
                X = HaltonPoints(:,1);
                Y = HaltonPoints(:,2);
                ind = X.^2+Y.^2<1-h/2;
                Xi = X(ind);
                Yi = Y(ind);
            case 'Uniform'
                N1 = ceil((b-a)/h);
                N2 = ceil((d-c)/h);
                yh=-1:2/N1:1;
                yv=-1:2/N2:1;
                [X, Y] = meshgrid(yh',yv');
                X = X(:);
                Y = Y(:);
                ind = X.^2+Y.^2<1-2/N1;
                Xi = X(ind);
                Yi = Y(ind);
        end
        Points{1,1} = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Points{1,2} = 'Interior points';
        Theta = (0:2/ceil((b-a)/h):2*pi)';
        Xi = ((b-a)/2).*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2).*sin(Theta)+(d+c)/2;
        Points{2,1} = [Xi,Yi];
        Points{2,2} = 'Boundary points';
    case 3 % Irregular Domain
        switch dist_type
            case 'Halton'
                N  = floor(((b-a)/h)*((d-c)/h));
                HaltonPoints = Halton(N,-1,1,-1,1);
                X = HaltonPoints(:,1);
                Y = HaltonPoints(:,2);
            case 'Uniform'
                N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
                yh=-1:2/N1:1;
                yv=-1:2/N2:1;
                [X, Y] = meshgrid(yh',yv');
        end
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =(0.7 + 0.12*(sin(6*Theta)+sin(3*Theta)));
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Points{1,1} = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Points{1,2} = 'Interior points';
        Theta = (0:h:2*pi)';
        R = (0.7 + 0.12*(sin(6*Theta)+sin(3*Theta)));
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Points{2,1} = [Xi,Yi];
        Points{2,2} = 'Boundary points';
    case 4 % Irregular Domain
        switch dist_type
            case 'Halton'
                N  = floor(((b-a)/h)*((d-c)/h));
                HaltonPoints = Halton(N,-1,1,-1,1);
                X = HaltonPoints(:,1);
                Y = HaltonPoints(:,2);
            case 'Uniform'
                N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
                yh=-1:2/N1:1;
                yv=-1:2/N2:1;
                [X, Y] = meshgrid(yh',yv');
        end
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =(.5 + .25*(sin(4*Theta).^2+sin(5*Theta)))+.2*sin(Theta).^3;
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Points{1,1} = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Points{1,2} = 'Interior points';
        Theta = (0:h:2*pi)';
        R = (.5 + .25*(sin(4*Theta).^2+sin(5*Theta)))+.2*sin(Theta).^3;
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Points{2,1} = [Xi,Yi];
        Points{2,2} = 'Boundary points';
    case 5 % Irregular Domain
        switch dist_type
            case 'Halton'
                N  = floor(((b-a)/h)*((d-c)/h));
                HaltonPoints = Halton(N,-1,1,-1,1);
                X = HaltonPoints(:,1);
                Y = HaltonPoints(:,2);
            case 'Uniform'
                N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
                yh=-1:2/N1:1;
                yv=-1:2/N2:1;
                [X, Y] = meshgrid(yh',yv');
        end
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =(.55 + .25*(sin(4*Theta).^2+sin(6*Theta)));
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Points{1,1} = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Points{1,2} = 'Interior points';
        Theta = (0:h:2*pi)';
        R = (.55 + .25*(sin(4*Theta).^2+sin(6*Theta)));
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Points{2,1} = [Xi,Yi];
        Points{2,2} = 'Boundary points';
    case 6 % Irregular Domain
        switch dist_type
            case 'Halton'
                N  = floor(((b-a)/h)*((d-c)/h));
                HaltonPoints = Halton(N,-1,1,-1,1);
                X = HaltonPoints(:,1);
                Y = HaltonPoints(:,2);
            case 'Uniform'
                N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
                yh=-1:2/N1:1;
                yv=-1:2/N2:1;
                [X, Y] = meshgrid(yh',yv');
        end
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =0.25*(2 + sin(2*Theta)-0.01*cos(5*Theta-pi/2)+0.63*sin(6*Theta-.1));
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Points{1,1} = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Points{1,2} = 'Interior points';
        Theta = (0:h:2*pi)';
        R = 0.25*(2 + sin(2*Theta)-0.01*cos(5*Theta-pi/2)+0.63*sin(6*Theta-.1));
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Points{2,1} = [Xi,Yi];
        Points{2,2} = 'Boundary points';
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



