function sigma = compute_surftension(rho,E,idx,eqstate)
%COMPUTE_SURFTENSION
%Compute surface tension coefficient sigma due to surfactants, either by
%linear or nonlinear equation of state
%
%   sigma = compute_surftension(rho,E,idx)
%
%Returns:
%  **sigma** -- surface tension coefficient
%
%:param rho: surfactant concentration (non-dimensionalized)
%:param E: elasticity parameter between 0 and 1
%

sigma = zeros(size(rho));
for j=1:size(idx,1)
    I = idx(j,1):idx(j,2);
    
    switch eqstate
        case 'linear'
            %Linear equation of state, only below dilute concentrations
            sigma(I) = 1-E(j)*rho(I);
        case 'langmuir'
            %Non-linear equation of state (Langmuir)
            sigma(I) = 1+E(j)*log(1-rho(I));
    end
end