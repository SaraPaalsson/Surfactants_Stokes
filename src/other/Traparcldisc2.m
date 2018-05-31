function [z,W] = Traparcldisc2(f,N,NL,pointsPerPanel)
%[z,W] = Traparcldisc(f,N)
%[z,W] = Traparcldisc(f,N,NL)
%
%Discretizes the curve given by f(t) t = [0,2*pi] using the trapezoidal rule
%equidistant in arc-length .
%N is the number of quadrature points requested and if specified,
%NL is the number of 16-point Gauss-Legendre panels used to compute the 
%length of the curve. Default is NL = 1000.
%Non-adaptive. Make sure that the curve is well resolved by N points.

%if nargin < 5
if NL == 0 %SARATEST %CHECK
    NL = 1000;
end
if pointsPerPanel == 16
    [T,W] = GLinterval(0,2*pi,NL);
else
    [T,W] = GLinterval32(0,2*pi,NL);
end
[~,zp] = f(T);
L = W'*abs(zp);

L2 = linspace(1/N,1-1/N,N-1)'*L;
%Initial guess
t = linspace(0,2*pi,N+1)';
t = t(2:end-1);

dt = 1;
iter = 0;
while norm(dt)/norm(t) > 1e-13 && iter < 30
    %Set up 16-point GL quadrature on each segment.
    if pointsPerPanel == 16
        [T,W] = GLinterval(0,2*pi,N-1,[0;t]);
    else
        [T,W] = GLinterval32(0,2*pi,N-1,[0;t]);
    end
    [~,zp] = f([t;T]);
    %Compute the cumulative sum of all the segment lengths
    if pointsPerPanel == 16
        F = cumsum(sum(reshape(W.*abs(zp(N:end)),16,N-1))');
    else
        F = cumsum(sum(reshape(W.*abs(zp(N:end)),32,N-1))');
    end
    dt = (F-L2)./abs(zp(1:N-1));
    %Sort the parameters just in case. 
    t = sort(t - dt);
    
    iter = iter + 1;
end
if iter == 30
    warning('Traparcldisc : Newton did not converge, bad discretization');
end
t = [0;t];

z = f(t);
W = 2*pi/N*ones(N,1);
