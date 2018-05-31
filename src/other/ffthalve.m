function zout = ffthalve(z,idx)

zout = zeros(length(z)/2,1);
for j = 1:size(idx,1)
    i1 = 2*(idx(j,1)-1)+1;
    i2 = 2*idx(j,2);
    N = i2-i1+1;
    c = fft(z(i1:i2));
    zout(idx(j,1):idx(j,2)) = 0.5*ifft([c(1:N/4);c(3*N/4+1:end)]);
end