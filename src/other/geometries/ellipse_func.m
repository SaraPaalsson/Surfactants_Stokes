function [z,zp,zpp] = ellipse_func(a,b,center,theta,phi)

% t = pi/2; %For crowdy valid!
t = 0;

s = sin(phi+t);
c = cos(phi+t);

%Calculate z,zp and zpp at the nodes.
ell = (a*c + 1i*b*s)*exp(1i*theta);
z = center + ell;
zp = (1i*b*c - a*s)*exp(1i*theta);
zpp = -ell;
