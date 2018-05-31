function [z,zp,zpp] = slantedflower(a,N,r,theta,t)
% zFuncs{1} = @(T) slantedflower(0.6,6,0.4,2,T);
z = exp(1i*(theta+t)).*(1+a*cos(N*t)).*(1+r*cos(t));
zp = 1i*exp(1i*(theta+t)).*(a*r*cos(t).*cos(N*t)+1i*a*N*r*cos(t).*sin(N*t)+...
    1i*a*r*sin(t).*cos(N*t) + 1i*a*N*sin(N*t)+a*cos(N*t) + 1i*r*sin(t)+r*cos(t) +1);
zpp = -exp(1i*(theta+t)).*(a*N^2*r*cos(t).*cos(N*t)+a*N^2*cos(N*t)-2*a*N*r*sin(t).*sin(N*t)+...
    2*a*r*cos(t).*cos(N*t) + 2i*a*N*r*cos(t).*sin(N*t)+2i*a*r*sin(t).*cos(N*t)+...
    2i*a*N*sin(N*t)+a*cos(N*t)+2i*r*sin(t)+2*r*cos(t)+1);