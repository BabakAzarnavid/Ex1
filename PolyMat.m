function p = PolyMat(X,op)
% This function computes the bisis polynomials of degree at most "m" in
% "d" dimensional space at "X" which "X" is the set of "n" points in this
% space.
%Example:  Let X={[1 2],[3 4],[5 6]} is set of two points in the plan then
% the basis polynomials of degree 2 in this space is:
%   f=[ 1, x, y, x^2, x*y, y^2],
%  so f(X)=
%       1     1     2     1     2     4
%       1     3     4     9    12    16
%       1     5     6    25    30    36
%
global Npoly
 Dpoly =Npoly-1;
if nargin ==1
    op ='0';
end
if Dpoly<0
p=[];
else
[n,d] = size(X);
q=nchoosek(Dpoly+d,d);
A=zeros(n,q);
k=1;
if strcmp( op,'0')==1
    for i=0:Dpoly
        if d==2
            for j=0:i
                A(:,k)=(X(:,1)).^(i-j).*(X(:,2)).^(j);
                k=k+1;
            end
        end
        p=A;
    end
else
    [h1,h2]=size(op);
    if h2>d
        disp('error, the number of variables are less than partial derivative order');
        return;
    end
    for i=0:Dpoly
        if d==1
            if strcmp(op,'x')==1
                if i<1
                    A(:,k)=0;
                else
                    A(:,k)=i*X.^(i-1);
                end
            elseif strcmp(op,'xx')==1
                if i<2
                    A(:,k)=0;
                else
                    A(:,k)=i*(i-1)*X.^(i-2);
                end
            elseif strcmp(op,'xxx')==1
                if i<3
                    A(:,k)=0;
                else
                    A(:,k)=i*(i-1)*(i-2)*X.^(i-3);
                end
            end
            k=k+1;
        end
        %A(:,k)=(X(:,1)).^(i-j).*(X(:,2)).^(j);
        if d==2
            for j=0:i
                switch op
                    case 'x'
                        if i-j<1
                            A(:,k)=0;
                        else
                            A(:,k)=(i-j)*(X(:,1)).^(i-j-1).*(X(:,2)).^(j);
                        end
                    case 'xx'
                        if i-j<2
                            A(:,k)=0;
                        else
                            A(:,k)=(i-j)*(i-j-1)*(X(:,1)).^(i-j-2).*(X(:,2)).^(j);
                        end
                    case 'xy'
                        if (i-j<1) || (j<1)
                            A(:,k)=0;
                        else
                            A(:,k)=j*(i-j)*(X(:,1)).^(i-j-1).*(X(:,2)).^(j-1);
                        end
                    case 'yx'
                        if (i-j<1) || (j<1)
                            A(:,k)=0;
                        else
                            A(:,k)=j*(i-j)*(X(:,1)).^(i-j-1).*(X(:,2)).^(j-1);
                        end
                        %A(:,k)=(X(:,1)).^(i-j).*(X(:,2)).^(j);
                    case'y'
                        if j<1
                            A(:,k)=0;
                        else
                            A(:,k)=(j)*(X(:,1)).^(i-j).*(X(:,2)).^(j-1);
                        end
                    case 'yy'
                        if j<2
                            A(:,k)=0;
                        else
                            A(:,k)=(j)*(j-1)*(X(:,1)).^(i-j).*(X(:,2)).^(j-2);
                        end
                    case 'L'
                        if i-j<2
                            A(:,k)=0;
                            Lx = A(:,k);
                        else
                            A(:,k)=(i-j)*(i-j-1)*(X(:,1)).^(i-j-2).*(X(:,2)).^(j);
                            Lx = A(:,k);
                        end
                        if j<2
                            A(:,k)=0;
                            Ly = A(:,k);
                        else
                            A(:,k)=(j)*(j-1)*(X(:,1)).^(i-j).*(X(:,2)).^(j-2);
                            Ly = A(:,k);
                        end
                        A(:,k) = Lx+Ly;
                end
                k=k+1;
            end
        end
        %A(:,k)=X(:,1).^(i-j-l).*X(:,2).^(l).*X(:,3).^j;
        if d==3
            for j=0:i
                for l=0:i-j
                    switch op
                        case 'x'
                            if i-j-l<1
                                A(:,k)=0;
                            else
                                A(:,k)=(i-j-l)*X(:,1).^(i-j-l-1).*X(:,2).^(l).*X(:,3).^j;
                            end
                        case 'y'
                            if l<1
                                A(:,k)=0;
                            else
                                A(:,k)=(l)*X(:,1).^(i-j-l).*X(:,2).^(l-1).*X(:,3).^j;
                            end
                        case 'z'
                            if j<1
                                A(:,k)=0;
                            else
                                A(:,k)=j*X(:,1).^(i-j-l).*X(:,2).^(l).*X(:,3).^(j-1);
                            end                            
                        case 'xx'
                            if i-j-l<2
                                A(:,k)=0;
                            else
                                A(:,k)=(i-j-l)*(i-j-l-1)*X(:,1).^(i-j-l-2).*X(:,2).^(l).*X(:,3).^j;
                            end
                        case 'yy'
                            if l<2
                                A(:,k)=0;
                            else
                                A(:,k)=(l)*(l-1)*X(:,1).^(i-j-l).*X(:,2).^(l-2).*X(:,3).^j;
                            end
                        case 'zz'
                            if j<2
                                A(:,k)=0;
                            else
                                A(:,k)=j*(j-1)*X(:,1).^(i-j-l).*X(:,2).^(l).*X(:,3).^(j-2);
                            end
                        case 'L'
                            if i-j-l<2
                                A(:,k)=0;
                                Lx = A(:,k);
                            else
                                A(:,k)=(i-j-l)*(i-j-l-1)*X(:,1).^(i-j-l-2).*X(:,2).^(l).*X(:,3).^j;
                                Lx = A(:,k);
                            end
                            if l<2
                                A(:,k)=0;
                                Ly = A(:,k);
                            else
                                A(:,k)=(l)*(l-1)*X(:,1).^(i-j-l).*X(:,2).^(l-2).*X(:,3).^j;
                                Ly = A(:,k);
                            end
                            if j<2
                                A(:,k)=0;
                                Lz = A(:,k);
                            else
                                A(:,k)=j*(j-1)*X(:,1).^(i-j-l).*X(:,2).^(l).*X(:,3).^(j-2);
                                Lx = A(:,k);
                            end
                            A(:,k) = Lx+Ly+Lz;
                        case 'xy'
                            if i-j-l<1 || l<1
                                A(:,k)=0;
                            else
                                A(:,k)=(l)*(i-j-l)*X(:,1).^(i-j-l-1).*X(:,2).^(l-1).*X(:,3).^j;
                            end
                        case 'yx'
                            if i-j-l<1 || l<1
                                A(:,k)=0;
                            else
                                A(:,k)=(l)*(i-j-l)*X(:,1).^(i-j-l-1).*X(:,2).^(l-1).*X(:,3).^j;
                            end
                        case 'yz'
                            if l<1 || j<1
                                A(:,k)=0;
                            else
                                A(:,k)=(l)*j*X(:,1).^(i-j-l).*X(:,2).^(l-1).*X(:,3).^(j-1);
                            end
                        case 'zy'
                            if l<1 || j<1
                                A(:,k)=0;
                            else
                                A(:,k)=(l)*j*X(:,1).^(i-j-l).*X(:,2).^(l-1).*X(:,3).^(j-1);
                            end
                        case 'xz'
                            if j<1 || i-j-l<1
                                A(:,k)=0;
                            else
                                A(:,k)=j*(i-j-l)*X(:,1).^(i-j-l-1).*X(:,2).^(l).*X(:,3).^(j-1);
                            end
                        case 'zx'
                            if j<1 || i-j-l<1
                                A(:,k)=0;
                            else
                                A(:,k)=j*(i-j-l)*X(:,1).^(i-j-l-1).*X(:,2).^(l).*X(:,3).^(j-1);
                            end
                    end
                    k=k+1;
                end
            end
        end
        p=A;
    end
end
end