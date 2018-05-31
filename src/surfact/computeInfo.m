function [sa,ka] = computeInfo(z,idx,W)
%COMPUTEINFO
%Compute arclength and curvature on equidistant grid
% 
%  [sa,ka] = computeInfo(z,idx,W)
%
%Returns:
%  **sa** -- equal arclength parameter
%
%  **ka** -- curvature
%
%:param z: interface discretization points
%:param idx: index vector for multiple droplets
%:param W: quadrature weights
%

[zp,zpp] = fft_diff(z,idx);
curL = zeros(size(idx,1),1); sa = zeros(size(z)); ka = sa;
for k=1:size(idx,1)
    curL(k) = sum(W(idx(k,1):idx(k,2)).*abs(zp(idx(k,1):idx(k,2))));
    sa(idx(k,1):idx(k,2)) = curL(k)/(2*pi);
    ka(idx(k,1):idx(k,2)) = imag(zpp(idx(k,1):idx(k,2))./zp(idx(k,1):idx(k,2))) ...
        ./abs(zp(idx(k,1):idx(k,2)));
end