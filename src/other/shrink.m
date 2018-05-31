function z = shrink(z,N)

M = length(z);
tau = 10/M^2;
Mr = 2*M;

T = linspace(0,2*pi,N+1)';
T = T(1:end-1);
k = (-M/2:M/2-1)';
Fm = sqrt(pi/tau)*exp(k.^2*tau).*fftshift(fft(z))/M;

Fmtau = zeros(Mr,1);
Fmtau(1:M/2) = Fm(M/2+1:M);
Fmtau(Mr-M/2+1:Mr) = Fm(1:M/2);
kk = zeros(Mr,1);
kk(1:M/2) = k(M/2+1:M);
kk(Mr-M/2+1:Mr) = k(1:M/2);

fmtau = ifft(Fmtau)*Mr;
tmp = kk.*Fmtau;
fmtaup = ifft(1i*tmp)*Mr;
fmtaupp = ifft(-kk.*tmp)*Mr;

if sum(imag(z)) == 0
    [z,zp,zpp] = mex_Trap2GL_real(fmtau,fmtaup,fmtaupp,T);
else
    [z,zp,zpp] = mex_Trap2GL(fmtau,fmtaup,fmtaupp,T);
end