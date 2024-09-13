function A = Radialfun(Xe,X,rbftype,delta,op)
% rbfpar for wendland functions is Kind of wendland function
% rbfscale is the scale of trial function
% operator is the order of partial derivative of function
%rbf_kind = 'Wendland','Gaussian', 'MQ', 'IMQ', 'Polyharmonic'
global rbfpar
 n = size(Xe,1);
 [m,dim] = size(X);
 r  = distance_matrix(Xe,X)+eps;
 x = diff_matrix(Xe,X,'x');
 y = diff_matrix(Xe,X,'y');
%%
switch rbftype
    case 'Wendland30' % phi_3,0 in C^0
        e = 1/delta;
        y1 = max(1-e*r,0);
        switch op
            case '0'
                A = y1.^2;
            case 'x'
                A = NaN;
            case 'y'
                A = NaN;
            case 'L'
                A = NaN;
        end
    case 'Wendland31' % phi_3,1  in C^2
        e = 1/delta;
        y1 = max(1-e*r,0);
        switch op
            case '0'
                A = y1.^4.*(4*e*r+1);
            case 'x'
                A = -20*e^2.*y1.^3.*x;
            case 'y'
                A = -20*e^2.*y1.^3.*y;
            case 'L'
                phi_rDivr = -20*e^2.*y1.^3;
                phi_rr = 20*e^2*y1.^2.*(4*e*r-1);
                A = phi_rr + phi_rDivr;
        end
    case 'Wendland32'% phi_3,2   in C^4
        e = 1/delta;
        y1 = max(1-e*r,0);
        switch op
            case '0'
                A = y1.^6.*(35*(e*r).^2+18*(e*r)+3);
            case 'x'                
                A = -56*e^2*y1.^5.*(5*e*r+1).*x;
            case 'y'
                A = -56*e^2*y1.^5.*(5*e*r+1).*y;
            case 'L'
                phi_rDivr = -56*e^2.*y1.^5.*(5*e*r+1);
                phi_rr = 56*e^2*y1.^4.*(35*(e*r).^2-4*e*r-1);
                A = phi_rr + phi_rDivr;
        end
    case 'Wendland33'% phi_3,3   in C^6
        e = 1/delta;
        y1 = max(1-e*r,0);
        switch op
            case '0'
                A = y1.^8.*(32*(e*r).^3+25*(e*r).^2+8*(e*r)+1);
            case 'x'                
                A = -22*e^2*y1.^7.*(16*(e*r).^2+7*e*r+1).*x;
            case 'y'
                A = -22*e^2*y1.^7.*(16*(e*r).^2+7*e*r+1).*y;
            case 'L'
                phi_rDivr =  -22*e^2*y1.^7.*(16*(e*r).^2+7*e*r+1);
                phi_rr =  22*e^2*y1.^6.*(160*(e*r).^3+15*(e*r).^2-6*e*r-1);
                A = phi_rr + phi_rDivr;
        end
    case 'Wendland50' % % phi_5,0  in C^0
        e = 1/delta;
        y1 = max(1-e*r,0);
        switch op
            case '0'
                A = y1.^3;
            case 'x'
                A = NaN;
            case 'y'
                A = NaN;
            case 'L'                
                A = NaN;
        end
    case 'Wendland51'% phi_5,1  in C^2
        e = 1/delta;
        y1 = max(1-e*r,0);
        switch op
            case '0'
                A = y1.^5.*(5*e*r+1);
            case 'x'                
                A = -30*e^2.*(y1.^4).*x;
            case 'y'
                A = -30*e^2.*(y1.^4).*y;
            case 'L'
                phi_rDivr =  -30*e^2.*(y1.^4);
                phi_rr =  30*e^2*y1.^3.*(5*(e*r)-1);
                A = phi_rr + phi_rDivr;
        end
    case 'Wendland52'% phi_5,2    in C^6
        e = 1/delta;
        y1 = max(1-e*r,0);        
        switch op
            case '0'
                A = y1.^7.*(16*(e*r).^2+7*e*r+1);
            case 'x'
                A = -24*e*y1.^6.*(6*e*r+1).*x;
            case 'y'
                A = -24*e*y1.^6.*(6*e*r+1).*y;
            case 'L'
                phi_rDivr = -24*e*y1.^6.*(6*e*r+1);
                phi_rr = (20/delta^2)*y1.^2.*(4*r/delta-1);
                A = phi_rr + phi_rDivr./r;                
        end
    case 'GIMQ' % Generalized IMQ RBF
         e = 1/delta;
        switch op           
            case '0'
                A = 1./(1+(e*r).^2).^2;
        end
    case 'Matern'
        s = shape_par;
        e = 1/delta;
        b = -1/delta;
        switch op
            case '0'
                r=r+eps;
                A = (e*r).^(s-3/2).*besselk(s-3/2,e*r);
            case 'L' % Lplace = -2t*phi'+(1-t^2)phi''
                r = r+eps;
                yt = e^2*(e*r).^(s-3/2-1).*besselk(s-3/2-1,e*r);
                ytt = e^4*(e*r).^(s-3/2-2).*besselk(s-3/2-2,e*r);
                A = (1-t.^2).*ytt-2*t.*yt;
        end
    case 'Polyharmonic'% Polyharmonic spline
        k = rbfpar;
        if mod(k,2)==0
            r = r + eps;
            switch op
                case '0'
                    A = r.^k.*log(r);
                case 'x'
                    A = x.*r.^(k-2).*(k*log(r)+1);
                case 'y'
                    A = y.*r.^(k-2).*(k*log(r)+1);
                case 'L'
                    A = k*r.^(k-2).*(k*log(r)+2);
            end
        else
            r = r + eps;
            switch op
                case '0'
                    A = r.^k;
                case 'x'
                    A = k*x.*r.^(k-2);
                case 'y'
                    A = k*y.*r.^(k-2);
                case 'L'
                    A = (k^2+k*(dim-2))*r.^(k-2);
                case 'L2'
                    A = (k^2+k*(dim-2))*((k-2)^2+(k-2)*(dim-2))*r.^(k-4);
            end
        end
    case 'Gaussian'
        e = delta;
        switch op
            case '0'
                A = exp(-(e*r).^2);
            case 'x'
                phi_rDivr = -2*e^2*exp(-(e*r).^2);                
                A = x.*phi_rDivr;
            case 'y'
                phi_rDivr = -2*e^2*exp(-(e*r).^2);
                A = y.*phi_rDivr;
            case 'L'
                phi_rDivr = -2*e^2*exp(-(e*r).^2);
                phi_rr = 2*e^2.*exp(-(e*r).^2).*(2*(e*r).^2-1);
                A= phi_rr + (dim-1)*phi_rDivr; %uncompleted
        end
    case 'MQ'
        e = delta;
        switch op
            case '0'
                A = sqrt(1+(e*r).^2);
            case 'x'
                phi_rDivr = (e^2)./sqrt(1+(e*r).^2);
                A = x.*phi_rDivr;
            case 'y'
                phi_rDivr = (e^2)./sqrt(1+(e*r).^2);
                A = y.*phi_rDivr;
            case 'L'
                phi_rDivr = (e^2)./sqrt(1+(e*r).^2);
                phi_rr = e^2./(1+(e*r).^2).^(3/2);
                A = phi_rDivr + phi_rr;
        end
    case 'IMQ'
        e = delta;
        switch op
            case '0'
                A = 1./sqrt(1+(e*r).^2);
            case 'x'                
                phi_rDivr = -e^2./sqrt(1+(e*r).^2).^(3/2);
                A = x.*phi_rDivr;
            case 'y'                
                phi_rDivr = -e^2./sqrt(1+(e*r).^2).^(3/2);
                A = y.*phi_rDivr;
            case 'L'
                phi_rDivr = -e^2./sqrt(1+(e*r).^2).^(3/2);
                phi_rr = e^2*(2*(e*r).^2-1)./(1+(e*r).^2).^(5/2);
                A = phi_rDivr + phi_rr;
        end
end


