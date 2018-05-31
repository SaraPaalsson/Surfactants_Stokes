function F = fft_primitive(f,idx)
F = zeros(size(f));
for j = 1:size(idx,1)
    N = idx(j,2)-idx(j,1)+1;
    c = fft(f(idx(j,1):idx(j,2)));
    k = [(0:N/2) (-N/2+1:-1)]';
    c = -1i*c./k;
    c(1) = 0;
    F(idx(j,1):idx(j,2)) = ifft(c);
end