function [A,L,cg] = get_area(zFuncs,N)
A = zeros(length(zFuncs),1);
L = zeros(length(zFuncs),1);
cg = zeros(length(zFuncs),1);
[T,W] = TrapInterval(-pi,pi,N);
for j = 1:length(zFuncs)
    [z,zp] = zFuncs{j}(T);
    %The area of the bubble.
    A(j) = sum(W.*real(z).*imag(zp));
    %The length of the bubble.
    L(j) = sum(W.*abs(zp));
    cg(j) = sum(W.*(real(z).^2/2+1i*real(z).*imag(z)).*imag(zp))/A(j);
end