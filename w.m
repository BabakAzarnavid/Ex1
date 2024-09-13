function y = w(X,Xc,rbftype,rcov,op)
K = size(Xc,1);
w0 = Radialfun(X,Xc,rbftype,rcov,'0');
wx = Radialfun(X,Xc,rbftype,rcov,'x');
wy = Radialfun(X,Xc,rbftype,rcov,'y');
wL = Radialfun(X,Xc,rbftype,rcov,'L');
sumw0 = repmat(sum(w0')',1,K);
sumwx = repmat(sum(wx')',1,K);
sumwy = repmat(sum(wy')',1,K);
sumwL = repmat(sum(wL')',1,K);
switch op
    
    case '0'
        y= w0./sumw0;
    case 'x'
        y = (wx.*sumw0-sumwx.*w0)./sumw0.^2;
    case 'y'
        y = (wy.*sumw0-sumwy.*w0)./sumw0.^2;
    case 'L'
        y = ((wL.*sumw0-w0.*sumwL).*sumw0-...
            2*(wx.*sumwx+wy.*sumwy).*sumw0+...
            2*w0.*(sumwx.^2+sumwy.^2))./sumw0.^3;
end


