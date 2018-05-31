function [Rmax,Rmin,D] = calc_deform(z,W,idx)
%CALC_DEFORM
%Computes maximum and minimum radius of the ellipse described by z, as
%well as the formation D=(Rmax-Rmin)/(Rmax+Rmin)
%
%OBS: This only works if the droplet has an alliptical shape. 
%
%  [Rmax,Rmin,D] = calc_deform(z,W,idx)
%
%Returns:
%  **Rmax** -- maximum radius of droplet
%
%  **Rmin** -- minimum radius of droplet
%
%  **D** -- deformation of droplet
%
%:param z: interface discretization points
%:param W: quadrature weights
%:param idx: index vector for multiple droplets
%

distanceZ = zeros(length(z),length(z));
for ii = 1:length(z)
    distanceZ(ii,:) = abs(z(ii)-z);
end
[rows,cols,~] = find(max(max(distanceZ))==distanceZ);

zp = fft_diff(z,idx);
curA = abs(sum(W.*real(z).*imag(zp)));
Rmax = abs(z(rows(1))-z(cols(1)))/2;
Rmin = curA/(pi*Rmax);

D = (Rmax-Rmin)/(Rmax+Rmin);