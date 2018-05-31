%%
tol = 1e-10;

% Testing  the surfactants 
% Read in set up data from setup file (1)
[globaldata,timedisc,drops,surfact,datastruct] = readInput('setup_test1',[]);

N = length(datastruct.W);

L = 2*pi; %exact length of circle radius 1
phi = linspace(0,2*pi,length(drops.z)+1)'; phi = phi(1:end-1);

% Test 1: arc length and curvature
sa = surfact.salpha; ka = surfact.kappa;
assert(abs(sa(1)-L/2/pi) <= tol,'Problem with computation of arc length: wrong value')
% assert(max(diff(sa)) <= tol,'Problem with computation of arc length: not constant')
assert(abs(-1-ka(1)) <= tol, 'Problem with computation of curvature: wrong value')
%Should we test something more advanced here? 

% Test 2: diffusion-less computations of surfactants, rho = 0
% sa = surfact.salpha; ka = surfact.kappa;
rho = zeros(N,1); 
U = ones(N,1); T = U; S = T; 
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,sa,ka,U,T,S,drops.idx);
assert(max(abs(rhoknew-0)) <= tol,'Problem with comp_surf: rho=0') 

% Test 3: diffusion-less computations of surfactants, stretching
% sa = surfact.salpha; ka = surfact.kappa;
S = zeros(N,1); T = S; 
U = sin(phi);
rho = cos(phi);
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,sa,ka,U,T,S,drops.idx);
tmp = 1/N*fftshift(fft(-(ka.*U.*rho)));
assert(max(abs(rhoknew-tmp)) <= tol,'Problem with comp_surf: stretching term') 

% Test 4: diffusion-less computations of surfactants, convection; T!=0
% sa = surfact.salpha; ka = surfact.kappa;
S = zeros(N,1); U = S; 
T = cos(phi);
rho = sin(phi);
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,sa,ka,U,T,S,drops.idx);
tmp = 1/N*fftshift(fft(-T.*cos(phi)./sa));
assert(max(abs(rhoknew-tmp)) <= tol,'Problem with comp_surf: convection term T') 

% Test 5: diffusion-less computations of surfactants, convection; S!=0
% sa = surfact.salpha; ka = surfact.kappa;
U = zeros(N,1); T = U; 
S = cos(phi);
rho = sin(phi);
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,sa,ka,U,T,S,drops.idx);
tmp = 1/N*fftshift(fft((cos(phi).^2-sin(phi).^2)./sa));
assert(max(abs(rhoknew-tmp)) <= tol,'Problem with comp_surf: convection term S') 

