function [ p ] = g( X,t )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[d1,d2]=size(X);
x=X(:,1);
y=X(:,2);
if d2~=2
    disp('error in dimension in f(x)');
    return;
end
% p = (3+2*t-2*(t.^2).*cos(x+y)).*sin(x+y)+2*t.*cos(x+y);
p = exp(t).*sin(x+y).*(2+exp(2*t).*sin(x+y).^2);
%exp(t)*sin(x+y)*(2+exp(2*t)*sin(x+y)^2)
% p = (2*t+1)*sin(x).*sin(y);
end

