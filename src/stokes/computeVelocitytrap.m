function [velocity,iter,nmodifs] = computeVelocitytrap(z,idx,bubble,D,beta,gamma,shear,extension)
% Simplified version of computeVelocity using Trapezoidal rule and no
% special quadrature.
%
% OBS. not updated.

zc = 3+3i; %DANGER

%Double the number of points. Helps with stability.
z = fftdouble(z,idx);
didx = zeros(size(idx));
for j = 1:size(idx,1)
    didx(j,1) = 2*(idx(j,1)-1)+1;
    didx(j,2) = 2*idx(j,2);
end
%Compute the derivatives, as well as outward using normal and angle
%derivative.
[zp,zpp] = fft_diff(z,didx);
W = 2*pi/length(z)*ones(length(z),1);
n = -1i*zp./abs(zp);

izppzp = imag(zpp./zp);

pe = [];
for j = 1:size(idx,1)
    pe = [pe;z(didx(j,1):16:didx(j,2));z(didx(j,1))];
end

%Interpolate to Gauss-Legendre points.
% [z,zp,zpp,W] = Trap2GL(z,idx);
curL = sum(W.*abs(zp));
%The multipoles require all points in the unit box [-0.5,0.5]x[-0.5,0.5].
%This is no real problem since all integral operators are scale invariant.
[z_sc,zp_sc] = scaletrans(z,zp,zpp);

%Compute special quadrature modifications, both for the integral equation
%if needed, and for the velocity evaluation.
% [modifs,imodifs,rmodifs,cidx] = mex_dospecquad2(z,zp,W,pe,bubble,idx);
% nmodifs = size(cidx,1);
nmodifs = 0;
cidx = 0;
imodifs = 0;
rmodifs = 0;
%The solution to the integral equation is known if beta == 0
if beta ~= 0
    rhs = -gamma/2*zp./abs(zp)-beta*shear/2*conj(z)-beta*(1-extension*1i*conj(z));
    %     rhs = -zp./abs(zp)/2;
    [omega,~,iter] = mygmresg_el(rhs,0,eps,200,0,@(x,varargin) mulfunc(x,z,zp,z_sc,zp_sc,zpp,W,beta,zc,cidx,imodifs));
else
    omega = -zp./abs(zp)/4;
    iter = 1;
end

%Use the density omega to compute the velocity of all points making up the
%interface.
omegaprime = fft_diff(omega,didx);
% figure
% T = TrapInterval(0,2*pi,length(z));
% y = mulfunc(rhs,z,zp,z_sc,zp_sc,zpp,W,beta,zc,cidx,imodifs);
% plot(T,real(omega),T,imag(omega))
% plot(T,real(y),T,imag(y))

% omegaprime = computederivative(omega,D,idx);
[~,M5] = mex_M1M5(z_sc,zp_sc,omega.*W/pi,1e-13,8,8);
% M1 = real(mex_M1(z_sc,zp_sc.*W/pi,1e-13,8,8));
M1 = 0;
velocity =-omegaprime.*W/pi-W/2/pi.*real(zpp./zp).*omega + M1.*omega+ M5-conj(omega).*W.*imag(zpp.*conj(zp))./conj(zp).^2/2i/pi+imag(shear*z)+extension*conj(z);
% velocity =-omegaprime.*W/pi-W/2/pi.*real(zpp./zp).*omega+ M5-conj(omega).*W.*imag(zpp.*conj(zp))./conj(zp).^2/2i/pi+imag(shear*z)+extension*conj(z);

% velocity =-omegaprime.*W/pi-M1r-1i*M1i +M1.*omega+ M5-conj(omega).*W.*imag(zpp.*conj(zp))./conj(zp).^2/2i/pi+imag(shear*z)+extension*conj(z);
% clf
% plot(T,real(omega),T,real(omegaprime))
% pause
% [modifs,nmodifs] = mex_dospecquad(z,zp,W,pe,bubble,idx,omega);
% velocity = velocity + modifs;

%Apply special quadrature corrections to the velocity
% velocity = velocity+omega.*rmodifs;
% mex_applymodifs(velocity,omega,modifs,cidx);
% velocity = velocity*curL/2/pi;
%Interpolate back to equidistant points
% velocity = GL2Trap(velocity);

%Maintain arc-length
U = real(velocity.*conj(n));
T = -real(fft_primitive(U.*izppzp,didx));
velocity = ffthalve((U+1i*T).*n,idx);

% velocity = ffthalve(velocity,idx);
end