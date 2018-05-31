function z = GL2Trap(z,pointsPerPanel)
%GL2TRAP
%Goes from Gauss-Legendre to equidistant points. Straight interpolation.
%
%  z = GL2Trap(z,pointsPerPanel)
%
%Returns:
%  **z** -- equidistant discretization points of interface
% 
%:param z: Gauss-Legendre nodes on panels
%:param pointsPerPanel: number of discretization points on each panel
%

if pointsPerPanel == 16
    load 'GL2TrapMat'
    for j = 1:length(z)/16
        z((j-1)*16+(1:16)) = I*z((j-1)*16+(1:16));
    end
else
%    for j = 1:length(z)/32
%       z((j-1)*16+(1:16)) = interp1(angle(z),z, 
%    end
end