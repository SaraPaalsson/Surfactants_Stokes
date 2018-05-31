function [velocity,iter,nmodifs,omega,U,T,S] =  ...
    computeVelocityOLD(z,idx,bubble,D,beta,gamma,shear,extension,omegaguess,sigma)
%COMPUTEVELOCITY
%Computes velocity through a boundary integral formulation of Stokes
%equations, both on droplet boundaries and in domain. Uses special
%quadrature to reduce errors and FMM to speed up computations.
%
%Returns: 
%   **velocity** --  the velocity with which to move the droplets
%
%   **iter** --  number of iterations needed for gmres
%
%   **nmodis** -- number of modifications by special quadrature
%
%   **omega** -- complex density computed by gmres
%
%   **U,T,S** -- normal and tangential velocities
%
%   **velocityDomain** -- velocity of domain points, [] if none
%
%:param z: boundary discretization
%:param idx: index table for multiple droplets
%:param bubble: bubble vector
%:param D: interpolation coefficients to go from G-L quadrature to equidistant
%:param beta,gamma: viscosity parameters
%:param shear,extension: imposed far-field velocities
%:param omegaguess: initial guess omega
%:param sigma: surface tension coefficient
%:param pointsPerPanel: number of discretization points per interface panel
%:param zDomain: domain discretization points
%

zc = zeros(size(idx,1),1);
for j=1:size(idx,1)
    zc(j) = sum(z(idx(j,1):idx(j,2)))/length(z(idx(j,1):idx(j,2)));
end



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
n = -1i*zp./abs(zp); %inward if new param
izppzp = imag(zpp./zp); %neg if new param

pe = [];
for j = 1:size(idx,1)
    pe = [pe;z(didx(j,1):16:didx(j,2));z(didx(j,1))];
end

% Save equidistant info for later velocity computations
zpeq = zp; %         

%Interpolate to Gauss-Legendre points.
[z,zp,zpp,W] = Trap2GL(z,didx,16);
zp = -zp; %Changed to new parametrisation 

%The multipoles require all points in the unit box [-0.5,0.5]x[-0.5,0.5].
%This is no real problem since all integral operators are scale invariant.
[z_sc,zp_sc] = scaletrans(z,zp,zpp);

%Compute special quadrature modifications, both for the integral equation
%if needed, and for the velocity evaluation.
zp_old = -zp+0; % Back to old parametrisation for spec quadrature
[modifs,imodifs,rmodifs,cidx] = mex_dospecquad3(z,zp_old,W,pe,bubble,idx,beta,4);
nmodifs = size(cidx,1);

sigmaD = fftdouble(sigma,idx);                                              %SURF%
[sigmaR,~,~,~] = Trap2GL_real(sigmaD,didx,16);                                   %SURF%
sigma = real(sigmaR);

%The solution to the integral equation is known if beta == 0
if norm(beta) ~= 0

    rhs = -gamma.*sigma.*zp./abs(zp)/2 + 1i*beta.*extension.*conj(z) - ...
        beta.*shear/2.*conj(z);
    
    [omega,~,iter] = mygmresg_el(rhs,0,1e-10,1000,0,@(x,varargin) mulfunc(x,z,zp,z_sc,zp_sc,zpp,W,beta,zc,cidx,imodifs,didx));

else
    omega = -sigma.*zp./abs(zp)/4;                                          %SURF%
    iter = 1;
end

%Use the density omega to compute the velocity of all points making up the
%interface.
omegaprime = -computederivative(omega,D,didx,16); %New paramaetrisation

[~,M5] = mex_M1M5(z_sc,zp_sc,omega.*W/pi,1e-13,8,8);
M1 = real(mex_M1(z_sc,zp_sc.*W/pi,1e-13,8,8));

velocity =-omegaprime.*W/pi + M1.*omega + M5-conj(omega).*W.*imag(zpp.*conj(zp))./conj(zp).^2/2i/pi;
velocity = velocity  + imag(shear*z)+extension*conj(z); %add far field


%Apply special quadrature corrections to the velocity
rmodifs = -rmodifs;
velocity = velocity+omega.*rmodifs;

modifs = -modifs;
mex_applymodifs(velocity,omega,modifs,cidx);


% %Interpolate back to equidistant points
velocity = GL2Trap(velocity,16);

for j=1:size(didx,1)
    I = didx(j,1):didx(j,2);
    N = length(z(I)); 
    k2 = (-N/2:N/2-1)'; %SARATEST
    vk = 1/N*fftshift(fft(velocity(I)));
    Ik = (abs(vk)>1e-12); %13? 12? 11? 10?
    vk = vk.*Ik;
    velocity(I) = N*ifft(ifftshift(vk));
end

%Maintain arc-length
S = real(velocity.*conj(zpeq./abs(zpeq)));                                  
U = real(velocity.*conj(n));
T = real(fft_primitive(U.*izppzp,didx)); % new parametrisation

velocity = ffthalve((U+1i*T).*n,idx);
U = ffthalve(U,idx);                                                        
T = ffthalve(T,idx);                                                        
S = ffthalve(S,idx);



end