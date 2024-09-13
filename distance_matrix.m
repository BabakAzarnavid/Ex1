function A = distance_matrix(X,Y)

% x=X*X';
% y=Y*Y';
% %normx = diag(x);
% %normy = diag(y);
% ux = ones(size(diag(x)));
% uy = ones(size(diag(y)));
% % xx = diag(x)*uy';
% % yy = ux*diag(y)';
% A = sqrt(diag(x)*uy'+ux*diag(y)'-2*X*Y');


%% Another algorithm
[n,dim]=size(X);
m = size(Y,1);
A = zeros(n,m);
for d=1:dim
A = A + (repmat(X(:,d),1,m)-repmat(Y(:,d),1,n)').^2;
end
A = sqrt(A);