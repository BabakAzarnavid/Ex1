function [Xall,Xin,Xbd] = points_sq_unity(a,b,c,d,h,type)
%UPOINTS Summary of this function goes here
%   Detailed explanation goes here
%
%   Ommega: the points in the region
%   Gamma1: the points on the inner boundary of the region
%   Gamma2= the points on the outer boundary
%   N= number of total points
switch type
    case 'R'
        yh=a-eps:h:b+eps;
        yv=c-eps:h:d+eps;
        
        %GL=[a*ones(n2+2,1),[c;yv';d]];
        
        
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
        %GR=[b*ones(n2+2,1),[c;yv';d]];
        GR=[b*ones(n2,1),yv'];
        GU=[yh(2:end-1)',d*ones(n1-2,1)];
        GD=[yh(2:end-1)',c*ones(n1-2,1)];
        Xbd=[GL;GR;GU;GD];
    case 'C'
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
        Xin = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Theta = (0:2/ceil((b-a)/h):2*pi)';
        Xi = ((b-a)/2).*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2).*sin(Theta)+(d+c)/2;
        Xbd = [Xi,Yi];
        Xall = [Xin;Xbd];
    case 'P3'
        N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
        yh=-1:2/N1:1;
        yv=-1:2/N2:1;
        [X, Y] = meshgrid(yh',yv');
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =(0.7 + 0.12*(sin(6*Theta)+sin(3*Theta)));
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Xin = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Theta = (0:h:2*pi)';
        R = (0.7 + 0.12*(sin(6*Theta)+sin(3*Theta)));
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Xbd = [Xi,Yi];
        Xall = [Xin;Xbd];
    case 'P4'
        N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
        yh=-1:2/N1:1;
        yv=-1:2/N2:1;
        [X, Y] = meshgrid(yh',yv');
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =(.5 + .25*(sin(4*Theta).^2+sin(5*Theta)))+.2*sin(Theta).^3;
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Xin = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Theta = (0:h:2*pi)';
        R = (.5 + .25*(sin(4*Theta).^2+sin(5*Theta)))+.2*sin(Theta).^3;
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Xbd = [Xi,Yi];
        Xall = [Xin;Xbd];
    case 'P5'
        N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
        yh=-1:2/N1:1;
        yv=-1:2/N2:1;
        [X, Y] = meshgrid(yh',yv');
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =(.55 + .25*(sin(4*Theta).^2+sin(6*Theta)));
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Xin = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2];
        Theta = (0:h:2*pi)';
        R = (.55 + .25*(sin(4*Theta).^2+sin(6*Theta)));
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Xbd = [Xi,Yi];
        Xall = [Xin;Xbd];
    case 'P6'
        N1 = ceil((b-a)/h);  N2 = ceil((d-c)/h);
        yh=-1:2/N1:1;
        yv=-1:2/N2:1;
        [X, Y] = meshgrid(yh',yv');
        X = X(:);
        Y = Y(:);
        [Theta,R] = cart2pol(X,Y);
        R0 =0.25*(2 + sin(2*Theta)-0.01*cos(5*Theta-pi/2)+0.63*sin(6*Theta-.1));
        ind = R<R0-h/8;
        Xi = X(ind);
        Yi = Y(ind);
        Xin = [((b-a)/2)*Xi+(b+a)/2,((d-c)/2)*Yi+(d+c)/2] ;
        Theta = (0:h:2*pi)';
        R = 0.25*(2 + sin(2*Theta)-0.01*cos(5*Theta-pi/2)+0.63*sin(6*Theta-.1));
        Xi = ((b-a)/2)*R.*cos(Theta)+(b+a)/2;
        Yi = ((d-c)/2)*R.*sin(Theta)+(d+c)/2;
        Xbd = [Xi,Yi];
        Xall = [Xin;Xbd];
        
end
