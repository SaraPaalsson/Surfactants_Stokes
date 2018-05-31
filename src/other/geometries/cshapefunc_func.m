function [z,zp,zpp] = cshapefunc_func(a,b,center,phi)

s = sin(phi);
c = cos(phi);
k = pi*(1-b);
e = -exp(-1i*c*k);

z = (s+a).*e + center;
zp = (c+k*1i*s.*(s+a)).*e;
zpp = -(s+k*(k*s.^2.*(s+a)-1i*c.*(a+3*s))).*e;



function [z,zp,zpp] = cshapefunc_funcold(a,b,center,phi)

s = sin(-phi);
c = cos(-phi);
e = exp(1i*c*pi*(1-b));

% z = (s+a).*e + center;
% zp = (c-(1-b)*pi*i*s.*(s+a)).*e;
% zpp = (-s -2*i*c.*s*pi*(1-b) - (s+a).*c*pi*(1-b)*i - (s+a).*s.^2*pi^2*(1-b)^2).*e;

z = (s+a).*e + center;
zp = -(c-(1-b)*pi*1i*s.*(s+a)).*e;
zpp = (-s - c*1i*pi*(1-b).*(3*s+a) - s.^2.*(s+a)*pi^2*(1-b)^2).*e;