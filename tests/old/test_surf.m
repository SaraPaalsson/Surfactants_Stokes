% Test surfactants

tol = 1e-10; %Absolute tolerance for tests

%Definitions 
N = 100; %Number disc. points
a = 1; b = 1; %Major, minor axes
c = 0; %Centre point
phi = linspace(0,2*pi,N+1)'; phi = phi(1:end-1);
[z,~,~] = ellipse_func(a,b,c,0,phi); 
W = 2*pi/N*ones(N,1);
idx = [1 N];

L = 2*pi; %exact length of circle radius 1

%% Test 1: arc length and curvature
[sa,ka] = computeInfo(z,idx,W);
assert(abs(sa(1)-L/2/pi) <= tol,'Problem with computation of arc length: wrong value')
assert(max(diff(sa)) <= tol,'Problem with computation of arc length: not constant')
assert(abs(ka(1)-1) <= tol, 'Problem with computation of curvature: wrong value')
%Should we test something more advanced here? 

%% Test 2: diffusion-less computations of surfactants, rho = 0
[sa,ka] = computeInfo(z,idx,W);
rho = zeros(N,1); 
U = ones(N,1); T = U; S = T; 
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,U,T,S,ka,sa,idx);
assert(max(abs(rhoknew-0)) <= tol,'Problem with comp_surf: rho=0') 

%% Test 3: diffusion-less computations of surfactants, stretching
[sa,ka] = computeInfo(z,idx,W);
S = zeros(N,1); T = S; 
U = sin(phi);
rho = cos(phi);
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,U,T,S,ka,sa,idx);
tmp = 1/N*fftshift(fft(-(ka.*U.*rho)));
assert(max(abs(rhoknew-tmp)) <= tol,'Problem with comp_surf: stretching term') 

%% Test 4: diffusion-less computations of surfactants, convection; T!=0
[sa,ka] = computeInfo(z,idx,W);
S = zeros(N,1); U = S; 
T = cos(phi);
rho = sin(phi);
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,U,T,S,ka,sa,idx);
tmp = 1/N*fftshift(fft(T.*cos(phi)./sa));
assert(max(abs(rhoknew-tmp)) <= tol,'Problem with comp_surf: convection term') 

%% Test 5: diffusion-less computations of surfactants, convection; S!=0
[sa,ka] = computeInfo(z,idx,W);
U = zeros(N,1); T = U; 
S = cos(phi);
rho = sin(phi);
rhok = 1/N*fftshift(fft(rho));
rhoknew = comp_surf(rhok,U,T,S,ka,sa,idx);
tmp = 1/N*fftshift(fft(-(cos(phi).^2-sin(phi).^2)./sa));
assert(max(abs(rhoknew-tmp)) <= tol,'Problem with comp_surf: convection term') 