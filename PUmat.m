function [A] = PUmat(op,X,Xe)
%P_UNITY Summary of this function goes here
%   Detailed explanation goes here
global  Xc rcov rbfw_type
K = size(Xc,1);
N = size(X,1);
M = size(Xe,1);
A = zeros(M,N);
W = w(Xe,Xc,rbfw_type,rcov,'0');
switch op
    case 'x'
        Wx = w(Xe,Xc,rbfw_type,rcov,'x');
    case 'y'
        Wy = w(Xe,Xc,rbfw_type,rcov,'y');
    case 'L'
        Wx = w(Xe,Xc,rbfw_type,rcov,'x');
        Wy = w(Xe,Xc,rbfw_type,rcov,'y');
        WL = w(Xe,Xc,rbfw_type,rcov,'L');
end
% NN = zeros(K,1);
type = 2; % 1 Classic , 2 Direct
switch type
    case 1%Classic
        for i=1:K
            ind_X    = find((Xc(i,1)-X(:,1)).^2 +  (Xc(i,2)-X(:,2)).^2<=rcov^2);
            ind_Xeval= find((Xc(i,1)-Xe(:,1)).^2+(Xc(i,2)-Xe(:,2)).^2<=rcov^2);
            if ~isempty(ind_Xeval)
                NN = length(ind_X);
                x  = X(ind_X,:);
                xe = Xe(ind_Xeval,:);
                xc = Xc(i,:);
                S = LagMat(xe,x,xc,0,'0');
                switch op
                    case '0'
                        % temp = repmat( W(ind_Xeval,i),1,NN).*S;
                        temp = repmat( W(ind_Xeval,i),1,NN).*S;
                        %A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                    case 'x'
                        Sx    = LagMat(xe,x,xc,1, 'x');
                        temp = repmat(W(ind_Xeval,i),1,NN).*Sx+repmat(Wx(ind_Xeval,i),1,NN).*S;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                    case 'y'
                        Sy    =  LagMat(xe,x,xc,1, 'y');
                        temp = repmat(W(ind_Xeval,i),1,NN).*Sy+repmat(Wy(ind_Xeval,i),1,NN).*S;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                    case 'L'
                        Sx    =  LagMat(xe,x,xc,1, 'x');
                        Sy    =  LagMat(xe,x,xc,1, 'y');
                        SL    =  LagMat(xe,x,xc,2, 'L');
                        temp1 = repmat(W(ind_Xeval,i),1,NN).*SL;
                        temp2 = 2*(repmat(Wx(ind_Xeval,i),1,NN).*Sx+repmat(Wy(ind_Xeval,i),1,NN).*Sy);
                        temp3 = repmat(WL(ind_Xeval,i),1,NN).*S;
                        %                temp = repmat( W(ind_Xeval,i),1,NN).*S;
                        %                  A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                        tempL = temp1+temp2+temp3;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+tempL;
                end
            end
        end
    case 2 %Direct
        for i=1:K
            ind_X    = find((Xc(i,1)-X(:,1)).^2 +  (Xc(i,2)-X(:,2)).^2<=rcov^2);
            ind_Xeval= find((Xc(i,1)-Xe(:,1)).^2+(Xc(i,2)-Xe(:,2)).^2<=rcov^2);
            if ~isempty(ind_Xeval)
                NN = length(ind_X);
                x  = X(ind_X,:);
                xe = Xe(ind_Xeval,:);
                xc = Xc(i,:);
                S = LagMat(xe,x,xc,0,'0');
                switch op
                    case '0'
                        % temp = repmat( W(ind_Xeval,i),1,NN).*S;
                        temp = repmat( W(ind_Xeval,i),1,NN).*S;
                        %A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                    case 'x'
                        Sx    = LagMat(xe,x,xc,1, 'x');
                        temp = repmat(W(ind_Xeval,i),1,NN).*Sx;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                    case 'y'
                        Sy    =  LagMat(xe,x,xc,1, 'y');
                        temp = repmat(W(ind_Xeval,i),1,NN).*Sy;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+temp;
                    case 'L'
                        SL    =  LagMat(xe,x,xc,2, 'L');
                        tempL = repmat(W(ind_Xeval,i),1,NN).*SL;
                        A(ind_Xeval,ind_X) = A(ind_Xeval,ind_X)+tempL;
                end
            end
        end
end
