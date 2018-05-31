function rhok = comp_surf(rhokold,sa,ka,U,T,S,idx)
%COMP_SURF
%Compute explicit part of surfactant ODE in Fourier space. 
%
%  rhok = comp_surf(rhokold,sa,ka,U,T,S,idx)
%
%Returns:
%  **rhok** -- fourier coefficients of function to be evaluated explicitly
%
%:param rhokold: fourier coeff. for previous timestep surfactant concentration
%:param sa: equal arclength parameter
%:param ka: curvature
%:param U,T,S: normal and tangential velocities
%:param idx: index vector for mutliple droplets
%

rhok = zeros(size(rhokold));
for j=1:size(idx,1)
    I = idx(j,1):idx(j,2);
    N = length(rhokold(I));

    %Zero pad
    M = 2*N;
    kM = (-M/2:M/2-1)';
    rhokM = [zeros((M-N)/2,1); rhokold(I); zeros((M-N)/2,1)];
    rhoM = real(ifft(M*ifftshift(rhokM)));
    Tk = fftshift(fft(T(I)));
    TkM = [zeros((M-N)/2,1); Tk; zeros((M-N)/2,1)];
    TM = ifft(M/N*ifftshift(TkM));
    Sk = fftshift(fft(S(I)));
    SkM = [zeros((M-N)/2,1); Sk; zeros((M-N)/2,1)];
    SM = ifft(M/N*ifftshift(SkM));
    Uk = fftshift(fft(U(I)));
    UkM = [zeros((M-N)/2,1); Uk; zeros((M-N)/2,1)];
    UM = ifft(M/N*ifftshift(UkM));
    kak = fftshift(fft(ka(I)));
    kakM = [zeros((M-N)/2,1); kak; zeros((M-N)/2,1)];
    kaM = ifft(M/N*ifftshift(kakM));

    %Differentiate rho
    drhokM = -1i*kM.*rhokM; %MINUS
    drhokM(abs(drhokM)<1e-12) = 0;
    drhoM = real(ifft(M*ifftshift(drhokM)));

    %Compute W1 
    W1 = TM.*drhoM./sa(j);
    w1k = fftshift(fft(W1));
    w1k(abs(w1k)<1e-12) = 0;

    %Compute W2
    Wtmp = rhoM.*SM./sa(j);
    wtmpk = fftshift(fft(Wtmp));
    wtmpk(abs(wtmpk)<1e-12) = 0;
    w2k = -1i*kM.*wtmpk; %MINUS

    %Compute W3
    W3 = kaM.*UM.*rhoM;
    w3k = fftshift(fft(W3));
    w3k(abs(w3k)<1e-12) = 0;

    rhoknewM = w1k-w2k-w3k;
    rhoknew = 1/M*rhoknewM((M-N)/2+1:(M+N)/2); %Downsample

    Ifilter = (abs(rhoknew)>1e-12); %Krasny filter
    frhok = rhoknew.*Ifilter;
    frhok(1) = 0; frhok(end) = 0;
    rhok(I) = frhok;
end

end
