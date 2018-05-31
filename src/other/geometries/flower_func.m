function [z,zp,zpp] = flower_func(a,scale,theta,Narms,center,phi)

z = scale*exp(1i*(phi+theta)).*(1+a*cos(Narms*phi)) + center;
zp = scale*exp(1i*(phi+theta)).*(1i + 1i*a*cos(Narms*phi) - a*Narms*sin(Narms*phi));
zpp = scale*exp(1i*(phi+theta)).*(-a*(1+Narms^2)*cos(Narms*phi)-1-2i*a*Narms*sin(Narms*phi));