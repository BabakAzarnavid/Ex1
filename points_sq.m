function [Ommega,Gamma,GL,GR,GU,GD] = points_sq( a,b,c,d,h)
yh=a+h:h:b-h;
yv=c+h:h:d-h;
n1=size(yh,2);
n2=size(yv,2);
%GL=[a*ones(n2+2,1),[c;yv';d]];
GL=[a*ones(n2,1),yv'];
%GR=[b*ones(n2+2,1),[c;yv';d]];
GR=[b*ones(n2,1),yv'];
GU=[yh',d*ones(n1,1)];
GD=[yh',c*ones(n1,1)];
Gamma=[GL;GR;GU;GD];
[Y1,Y2]= meshgrid(yh',yv');
Y=[Y1(:), Y2(:)];
Ommega=[Y(:,1),Y(:,2)];
end