%clc
tic;
clearvars
close all
global rbfscale rbf_type Npoly rbfpar do_scaling
global rbfw_type scaling_size
global X Xc Xeval h_y rcov
global h
a=0;b=1;c=0;d=1;
 alp=1/2;bt=1/2;
%alp=85/100;bt=4/10;
t_F= 1;
dt1 = t_F/10;
    niT1 = ceil(t_F/dt1);
    dist_type = 'Halton'; % Halton  Uniform
    domain_type = 1;
h=(b-a)/10;
nj=5;
    %% partition of unity options
    trial_func = 'Power';%  Power Wendland32  Wendland31 ...
    switch trial_func
        case 'Power'
            rbf_type = 'Polyharmonic';
            rbfpar = 10;
            Npoly = ceil(rbfpar/2);
            do_scaling = 1;
        case 'Wendland31' % Wendland32, Wendland33, Wendland32, Wendland51, Wendland52
            rbf_type = 'Wendland31';
            rbfscale = 4/(h);
            do_scaling = 0;
            rbfpar = 0;
            Npoly = ceil(rbfpar/2);
        case 'Gaussian'
            rbf_type = 'Gaussian';
            rbfscale = 1;
            do_scaling = 0;
            Npoly =0;
    end
    rbfw_type = 'Wendland32';
    h_y = 4*h;
     rcov =6*h_y ;
% rcov =sqrt(2)*h_y ;
        scaling_size = rcov;

    [Xc,Xci,Xcb] = points_sq_unity(a,b,c,d,h_y,'R');
    % [Xi,Xb,GL,GR,GU,GD]=points_sq(a,b,c,d,h);
    Points = Omega(a,b,c,d,h,dist_type,domain_type);
    Xi = Points{1,1};
    Xb = Points{2,1};
    X=[Xi;Xb];
    n=size(X,1);
    
    
    ni=length(Xi(:,1));
%     uui=zeros(ni,niT+1);
%     [xT,yT] = meshgrid(a:.3:b,c:.3:d);
%     Xeval = [xT(:) yT(:)];
    Points1 = Omega(a,b,c,d,h*sqrt(4),dist_type,domain_type);
    Xei = Points1{1,1};
    Xeb = Points1{2,1};
    Xeval=[Xei;Xeb];
    Aeval = PUmat('0',X,Xeval);
    Ai= PUmat('0',X,Xi);
    AL=PUmat('L',X,X);
    delAi= PUmat('L',X,Xi);
    %AG= PUmat('x',X,Xi)+PUmat('y',X,Xi);
    Ab = PUmat('0',X,Xb);

    

        %%%%article: 
% mu0=1;
% ex_sol=@(x,y,t) t.^(1+alp)*sin(pi*x).*sin(pi*y);
%  f0= @(x,y,t) t.^(1 + alp).*sin(pi*x).*sin(pi*y) + (1/gamma(2 + 2*alp)).*(2 *pi^2* t.^(1 + 2 *alp).* gamma(2 + alp).*sin(pi* x).*sin(pi* y)) + (1/gamma(2 + 2* alp + bt)).*(2*pi^2.*t.*(1 + 2 *alp + bt).* gamma(2 + alp).*sin(pi*x).*sin(pi*y)) + (1/gamma(3 + 3* alp)).*(t.^(2 + 3*alp).*gamma(3 + 2* alp).*(sin(pi*x).*sin(pi* y)).^2);
% g0=@(u) -u.^2;
% 
%  mu0=1;
%  ex_sol=@(x,y,t) t.^(1+alp)*sin(2*pi*x).*sin(2*pi*y);
% %%%%%alp=1/2;bt=1/2;
%  f0= @(x,y,t)  t.^(3/2).*sin (2*pi*x).*sin (2*pi*y) + 3* pi.^(5/2)* t.^2.*sin (2*pi*x).*sin (2*pi*y) +  16/5*pi^2*t.^(5/2).*sin (2*pi*x).*sin (2*pi*y) -  1/sqrt (pi)*(2*sqrt (t).* hypergeom ([1/3, 2/3, 1], [1/2, 1/2, 5/6, 7/6],  1/4*t.^3.*(sin (2*pi*x).*sin (2*pi*y)).^2) +     3/8*pi*t.^2.*      hypergeom ([5/6, 7/6], [1, 4/3, 5/3],  1/4*t.^3.*(sin (2*pi*x).*sin (2*pi*y)).^2).* sin (2*pi*x).*     sin (2*pi*y));
%       g0=@(u) exp(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 mu0=1;
 ex_sol=@(x,y,t) t.^(1+alp)*sin(2*pi*x).*sin(2*pi*y);
%%%%%alp=1/2;bt=1/2;
 f0= @(x,y,t)  t.^(3/2).*sin (2*pi*x).*sin (2*pi*y) +   3*pi^(5/2).*t.^2.*sin (2*pi*x).*sin (2*pi*y) +   16/5*pi^2.*t.^(5/2).* sin (2*pi*x).*sin (2*pi*y) -   3/8*sqrt (pi).*t.^2.*   hypergeom ([5/6, 7/6], [1, 4/3, 5/3], -1/4*       t.^3.* (sin (2*pi*x).*sin (2*pi*y)).^2) .*sin (2*pi*x).*   sin (2*pi*y);

%%%%%alp=85/100;bt=4/10;
 %f0= @(x,y,t) t.^(3/2).*sin (2*pi*x).*  sin (2*pi*y) + (1/(77 *gamma (5/4).*(gamma (7/4)).^2))*(18 *    sqrt (2).*pi.^(7/2).* (1./t).^(3/4).* t.^(7/2) .*sin (2*pi*x).*    sin (2*pi*y)) +  1/gamma (67/20)*(6*pi^(5/2).* t.^(47/20).*sin (2*pi*x).*    sin (2*pi*y)) -  1/(6* 3^(17/20).*gamma (67/60)* gamma (29/20).*     gamma (107/60))*(pi^(3/2).* t.^(47/20).*    hypergeom ([5/6, 7/6], [67/20, 29/20, 107/60], -1/4*        t.^3.*(sin (2*pi*x).*sin (2*pi*y)).^2) .*sin (2*pi*x).*    sin (2*pi*y));

 g0=@(u) sin(u);

%%%%%%%%%%%%%%%%%%%%%55
%  mu0=1;
%  ex_sol=@(x,y,t) t.^(1+alp)*sin(2*pi*x).*sin(2*pi*y);
% %%%%%alp=1/2;bt=1/2;
%  %f0= @(x,y,t)  t.^(3/2).*sin (2*pi*x).*sin (2*pi*y) +   3*pi^(5/2).*t.^2.*sin (2*pi*x).*sin (2*pi*y) +   16/5*pi^2.*t.^(5/2).* sin (2*pi*x).*sin (2*pi*y) -   3/8*sqrt (pi).*t.^2.*   hypergeom ([5/6, 7/6], [1, 4/3, 5/3], -1/4*       t.^3.* (sin (2*pi*x).*sin (2*pi*y)).^2) .*sin (2*pi*x).*   sin (2*pi*y);
% 
% %%%%%alp=1/2;bt=1/2;
%  f0= @(x,y,t) -1/sqrt(pi)*(2*sqrt(t).* hypergeom([1/3, 2/3, 1], [1/2, 1/2, 5/6, 7/6], -1/4*t.^3.*(sin(2*pi*x).*sin (2*pi*y)).^2)) + t.^(3/2).*sin (2*pi*x).*sin (2*pi*y) +  3*pi^(5/2)*t.^2.*sin (2*pi*x).*sin (2*pi*y) +  16/5*pi^2* t.^(5/2).*sin (2*pi*x).*sin (2*pi*y);
%  g0=@(u) cos(u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mu0=1;
%  ex_sol=@(x,y,t) t.^(1+alp)*sin(2*pi*x).*sin(2*pi*y);
% %%%%%alp=1/2;bt=1/2;
%  f0= @(x,y,t) -1/sqrt (pi)*(2 *sqrt (t).* hypergeom ([1/3, 2/3, 1, 1], [1/2, 5/6, 7/6], -t.^3.*(sin (2*pi*x) .*sin (2*pi*y)).^2)) + t.^(3/2).*sin (2*pi*x) .*sin (2*pi*y) + 3*pi^(5/2)*t.^2.*sin (2*pi*x) .*sin (2*pi*y) + 16/5*pi^2*t.^(5/2).*sin (2*pi*x) .*sin (2*pi*y);
%  g0=@(u) 1./(1+u.^2);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U0 = ex_sol(X(:,1),X(:,2),0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uu=zeros(n,niT1+1);

for k=1:1
    dt=dt1/k;
    niT = ceil(t_F/dt);
uu=zeros(n,niT+1);
er=zeros(n,niT+1);
 

uu(:,1)=U0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   COEFFICIENTS
syms x7
% omg=coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3)^(-alp),x7,'Order',niT+1),x7);
 omg=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3)^(-alp),x7,'Order',niT+1),x7));
% omg=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2)^(-alp),x7,'Order',niT+1),x7),300);
% omg=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3+(1-x7)^4/4)^(-alp),x7,'Order',niT+1),x7));
% omg=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3+(1-x7)^4/4+(1-x7)^5/5+(1-x7)^6/6)^(-alp),x7,'Order',niT+1),x7));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% omg2=coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3)^(-alp-bt),x7,'Order',niT+1),x7);

 omg2=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3)^(-alp-bt),x7,'Order',niT+1),x7));
% omg2=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2)^(-alp-bt),x7,'Order',niT+1),x7),300);
% omg2=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3++(1-x7)^4/4)^(-alp-bt),x7,'Order',niT+1),x7));
% omg2=vpa(coeffs(taylor(((1-x7)+(1-x7)^2/2+(1-x7)^3/3+(1-x7)^4/4+(1-x7)^5/5+(1-x7)^6/6)^(-alp-bt),x7,'Order',niT+1),x7));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%
i=1;

w11a=gamma(2)/gamma(alp+2)-omg(1);
w11b=gamma(2)/gamma(alp+bt+2)-omg2(1);


   B=zeros(n,n);
B(1:ni,:)=Ai-(mu0*dt^alp*(omg(1)+w11a)+dt^(alp+bt)*(omg2(1)+w11b))*delAi;
        B(ni+1:n,:) = Ab;
   %     rs=size(B);B5=B(:);bs=find(abs(B5)<10^(-8));B5(bs)=zeros(length(bs),1);B=reshape(B5,rs);
                        [L, U, P] = lu(B);
        Inv_Mass = U\(L\P)*eye(n);
 %      Inv_Mass = B\speye(size(B));
    

bh=mu0*dt^alp*sm(i,uu,omg,AL,alp)+dt^(alp+bt)*sm(i,uu,omg2,AL,alp+bt)+f0(X(:,1),X(:,2),i*dt)+dt^alp*smg(i,g0(uu),omg,alp);

for jj=1:nj    
   Rhs =bh+dt^alp*(omg(1)+w11a)*g0(uu(:,i+1));
   Rhs(ni+1:n) = zeros(size(Xb(:,1)));
   
%         uu(:,i) =B\[RhsI;Rhsb] ;
%  uu(:,i+1) =IB*[RhsI;Rhsb] ;
%  uu(:,i+1) =B\[RhsI;Rhsb] ;
 uu(:,i+1) =Inv_Mass*Rhs;
     %   uui=uu(1:ni,:);
          disp(i)
%       max(abs(uu(:,i+1)-ex_sol(X(:,1),X(:,2),i*dt)))
end
    er(:,i+1)=abs(uu(:,i+1)-ex_sol(X(:,1),X(:,2),i*dt));
    max(er(:,i+1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=2;

%%%%%%%%%%%
b07a=zeros(3,1);b07b=zeros(3,1);
%ac=[1 1 1 ; 0 1 2 ; 0 1 4];
aci=[1 -3/2 1/2; 0 2 -1; 0 -1/2 1/2];omg7a=omg(1:i+1);omg7b=omg2(1:i+1);

b07a(1)=1/gamma(alp+1)*i^alp-omg7a*ones(length(omg7a),1);
b07a(2)=gamma(2)/gamma(alp+2)*i^(1+alp)-flip(omg7a)*(0:length(omg7a)-1)';
b07a(3)=gamma(3)/gamma(alp+3)*i^(2+alp)-flip(omg7a)*((0:length(omg7a)-1).^2)';
w02a=aci*b07a;

b07b(1)=1/gamma(alp+bt+1)*i^(alp+bt)-omg7b*ones(length(omg7b),1);
b07b(2)=gamma(2)/gamma(alp+bt+2)*i^(1+alp+bt)-flip(omg7b)*(0:length(omg7b)-1)';
b07b(3)=gamma(3)/gamma(alp+bt+3)*i^(2+alp+bt)-flip(omg7b)*((0:length(omg7b)-1).^2)';
w02b=aci*b07b;
%%%%%%%%%%%


        B=zeros(n,n);
B(1:ni,:)=Ai-(mu0*dt^alp*(omg(1)+w02a(3))+dt^(alp+bt)*(omg2(1)+w02b(3)))*delAi;
        B(ni+1:n,:) = Ab;
   %     rs=size(B);B5=B(:);bs=find(abs(B5)<10^(-8));B5(bs)=zeros(length(bs),1);B=reshape(B5,rs);
                        [L, U, P] = lu(B);
        Inv_Mass = U\(L\P)*eye(n);
 %      Inv_Mass = B\speye(size(B));


bh=mu0*dt^alp*sm(i,uu,omg,AL,alp)+dt^(alp+bt)*sm(i,uu,omg2,AL,alp+bt)+dt^alp*smg(i,g0(uu),omg,alp)+f0(X(:,1),X(:,2),i*dt);

for jj=1:nj    
   Rhs =bh+dt^alp*(omg(1)+w02a(3))*g0(uu(:,i+1)); 
   Rhs(ni+1:n) = zeros(size(Xb(:,1)));
   
%         uu(:,i) =B\[RhsI;Rhsb] ;
%  uu(:,i+1) =IB*[RhsI;Rhsb] ;
%  uu(:,i+1) =B\[RhsI;Rhsb] ;
 uu(:,i+1) =Inv_Mass*Rhs ;
     %   uui=uu(1:ni,:);
     disp(i)

end
    er(:,i+1)=abs(uu(:,i+1)-ex_sol(X(:,1),X(:,2),i*dt));
max(er(:,i+1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% B1=Ai-(mu0*dt^alp*omg(1)+dt^(alp+bt)*omg2(1))*delAi;
%         B2 = Ab;
%         B = [B1;B2];
B=zeros(n,n);
B(1:ni,:)=Ai-(mu0*dt^alp*omg(1)+dt^(alp+bt)*omg2(1))*delAi;
        B(ni+1:n,:) = Ab;
   %     rs=size(B);B5=B(:);bs=find(abs(B5)<10^(-8));B5(bs)=zeros(length(bs),1);B=reshape(B5,rs);
                        [L, U, P] = lu(B);
        Inv_Mass = U\(L\P)*eye(n);
     %  Inv_Mass = B\speye(size(B));
     % Inv_Mass = sparse(B)\sparse(eye(size(B)));
%Inv_Mass =inv(B);
    for i = 3:niT
        
 
bh=mu0*dt^alp*sm(i,uu,omg,AL,alp)+dt^(alp+bt)*sm(i,uu,omg2,AL,alp+bt)+dt^alp*smg(i,g0(uu),omg,alp)+f0(X(:,1),X(:,2),i*dt);

for jj=1:nj    
   Rhs =bh+dt^alp*omg(1)*g0(uu(:,i+1));
   Rhs(ni+1:n) = zeros(size(Xb(:,1)));
   
%         uu(:,i) =B\[RhsI;Rhsb] ;
%  uu(:,i+1) =IB*[RhsI;Rhsb] ;
%  uu(:,i+1) =B\[RhsI;Rhsb] ;
 uu(:,i+1) =Inv_Mass*Rhs ;
     %   uui=uu(1:ni,:);
           disp(i)
 
end
er(:,i+1)=abs(uu(:,i+1)-ex_sol(X(:,1),X(:,2),i*dt));
max(er(:,i+1))
    end
    
end  
    
   
   err= max(er(:,i+1))
      erL2=sqrt(h^2*sum((er(:,i+1)).^2))
    %%
    err0= max(max(er))
      rel_err = norm(uu(:,i+1)-ex_sol(X(:,1),X(:,2),i*dt))/(norm(ex_sol(X(:,1),X(:,2),i*dt))+1)
% er1 = norm(uu2(:,i+1)-uu1(:,i+1))/(norm(uu2(:,i+1)));
    
%     ERR = abs(uu2(:,end)-uu1(:,end));
%    er2= max(ERR)
    % max(abs(uu(:,i+1)))
     
%     uapp=uu(:,i);
%clear Ab Ai delAi U1 V1 W1
%h=h/sqrt(2);


% erL2=sqrt(h^2*sum(ERR.^2))
% max(ERR)
%sqrt(h^2*sum(ERR.^2))
toc;
% 
%  S = scatteredInterpolant(X(:,1),X(:,2),ERR);sn0=20;
%  
%  figure('DefaultAxesFontSize',16); 
% % scatter3(X(:,1),X(:,2),S(X(:,1),X(:,2)),30,S(X(:,1),X(:,2)))
%  scatter(X(:,1),X(:,2),10,S(X(:,1),X(:,2)))
% % scatter(xb,yb,30,erb)
% 
% % scatter(xx3,yy3,30,S2(xx3,yy3))
% %scatter(xx3,yy3,30,S3(xx3,yy3))
% 
%  axis equal
% %axis([-1.6 2.8  -1.5  2])
% %axis([-.1 1.1   -.1 1.1])
% axis([0 1   0 1])
% %title('Approximate solution')
% title('Absolute error')
% 
% %grid on
% %grid minor
% xlabel('x')
% ylabel('y')
% 
% colormap('jet')
% colorbar 
% 
% 
% 
