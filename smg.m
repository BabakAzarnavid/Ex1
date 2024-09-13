function smm = smg(n1,u1,omga,alpha)
%%%% sm is approximation of --> (1/dt^(alpha))*I^(alpha) L1(u1) at t=n1*dt
%%%% using coffitients omega(bdf Lubish)

omga0=omga(1:n1+1);


if n1==1
    w0=zeros(3,1);
w0(1)=1/gamma(alpha+1)-gamma(2)/gamma(alpha+2)-omga(2);
w0(2)=gamma(2)/gamma(alpha+2)-omga(1);
w0(3)=0;

%%%
%%%%%%%%%%%%

smm=omga0(2)*u1(:,1)+w0(1)*u1(:,1);

%%%%%%%%%%%%

elseif n1==2

aci=[1 -3/2 1/2; 0 2 -1; 0 -1/2 1/2];
b00=zeros(3,1);


b00(1)=1/gamma(alpha+1)*n1^alpha-omga0*ones(length(omga0),1);
b00(2)=gamma(2)/gamma(alpha+2)*n1^(1+alpha)-flip(omga0)*(0:length(omga0)-1)';
b00(3)=gamma(3)/gamma(alpha+3)*n1^(2+alpha)-flip(omga0)*((0:length(omga0)-1).^2)';
w0=aci*b00;
%%%
%%%%%%%%%%%%

smm=(sum(((u1(:,1:2))*diag(flip(omga0(2:3))))'))'+(sum(((u1(:,1:2))*diag(w0(1:2)))'))';
%%%%%%%%%%%%

else

aci=[1 -3/2 1/2; 0 2 -1; 0 -1/2 1/2];
b00=zeros(3,1);
% i6=3;%%%%i>=2
b00(1)=1/gamma(alpha+1)*n1^alpha-omga0*ones(length(omga0),1);
b00(2)=gamma(2)/gamma(alpha+2)*n1^(1+alpha)-flip(omga0)*(0:length(omga0)-1)';
b00(3)=gamma(3)/gamma(alpha+3)*n1^(2+alpha)-flip(omga0)*((0:length(omga0)-1).^2)';
w0=aci*b00;
%%%%%%%%%%%%

smm=(sum(((u1(:,1:n1))*diag(flip(omga0(2:n1+1))))'))'+(sum(((u1(:,1:3))*diag(w0))'))';

%%%%%%%%%%%%
end



