function zout = fftdouble(z,idx)
zout = [];
for j = 1:size(idx,1)
    N = idx(j,2)-idx(j,1)+1;
    c = fft(z(idx(j,1):idx(j,2)));
    zout = [zout;2*ifft([c(1:N/2);zeros(N,1);c(N/2+1:end)])];
end