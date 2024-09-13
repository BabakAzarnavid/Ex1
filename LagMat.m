function L = LagMat(xe,x,xc,scaling_order, op)
global rbf_type rbfscale do_scaling scaling_size

if do_scaling
    xe = (xe-xc)/scaling_size;
    x = (x-xc)/scaling_size;
else
    xe = (xe-xc);
    x = (x-xc);
    scaling_order = 0;
end
A = Radialfun(x,x,rbf_type,rbfscale,'0');
P = PolyMat(x,'0');
q=size(P,2);
n = size(x,1);
Ae = Radialfun(xe,x,rbf_type,rbfscale,op)';
Pe = PolyMat(xe,op)';
L = [A P;P' zeros(q)]\[Ae;Pe];
L = (L(1:n,:)/scaling_size^scaling_order)';
end

