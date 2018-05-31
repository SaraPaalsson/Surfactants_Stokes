% Compute what initial surfactant concentration we can start with, given
% that the interface shrinks and we use the Langmuir eq.state
% 
% Created: 2018-04-18

close all; clear all; clc

% Interface at final time, i.e. when it is at its smallest
load('swissroll2_surf_outdata_t50.mat')
drops50 = simdata.drops;
surfact50 = simdata.surfact;

z50_1 = drops50.z(drops50.idx(1,1):drops50.idx(1,2));
W50_1 = drops50.W(drops50.idx(1,1):drops50.idx(1,2));

% What we want our surfactant concentration (uniform) to be at this small
% interface OBS this does not account for peaks!
M = 0.9;

rho50_1 = M*ones(size(z50_1));

% Compute area at final time
zp = fft_diff(z50_1,[1 length(z50_1)]);
curA = abs(sum(W50_1.*real(z50_1).*imag(zp)));
curL = sum(W50_1.*abs(zp));
sa50_1 = surfact50.salpha(1);

% The mass of surfactant at final time. This will be the same at initial
% time
rhomass_max = (sum(W50_1.*rho50_1*sa50_1));

% Interface at initial time, when it is at its largest
load('swissroll2_surf_outdata_t0.mat')
drops0 = simdata.drops;
surfact0 = simdata.surfact;

z0_1 = drops0.z(drops0.idx(1,1):drops0.idx(1,2));
W0_1 = drops0.W(drops0.idx(1,1):drops0.idx(1,2));
sa0_1 = surfact0.salpha(1);

% Compute area at time 50
zp = fft_diff(z0_1,[1 length(z0_1)]);
curA0 = abs(sum(W0_1.*real(z0_1).*imag(zp)));
curL0 = sum(W0_1.*abs(zp));

% Surfactant concentration at start
rho0 = rhomass_max/sum(sa0_1*W0_1)