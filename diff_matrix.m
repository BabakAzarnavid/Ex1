function A = diff_matrix(X,Y,op)
switch op
    case 'x'
        k=1;
    case 'y'
        k=2;
    case 'z'
        k=3;
    otherwise
        error('op is not valid')
end
A = repmat(X(:,k),1,size(Y,1))-repmat(Y(:,k),1,size(X,1))';