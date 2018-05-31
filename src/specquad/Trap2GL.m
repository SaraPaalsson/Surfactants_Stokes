function [z,zp,zpp,W] = Trap2GL(z,idx,tmp)
zp = zeros(size(z));
zpp = zeros(size(z));
W = zeros(size(z));
for j = 1:size(idx,1)
    i1 = idx(j,1);
    i2 = idx(j,2);
    
    M = i2-i1+1;
    tau = 10/M^2;
    Mr = 2*M;
    
    [T,W(i1:i2)] = GLinterval(0,2*pi,M/16);
    k = (-M/2:M/2-1)';
    Fm = sqrt(pi/tau)*exp(k.^2*tau).*fftshift(fft(z(i1:i2)))/M;
    
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
    
    [z(i1:i2),zp(i1:i2),zpp(i1:i2)] = mex_Trap2GL(fmtau,fmtaup,fmtaupp,T);
end