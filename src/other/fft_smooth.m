function z = fft_smooth(z,idx,ord)
%A Fourier filter of order ord
for j = 1:size(idx,1)
    N = idx(j,2)-idx(j,1)+1;
    coef = bandstop(N,ord,0.1);
    
    c = fft(z(idx(j,1):idx(j,2)));
    %     k = [(0:N/2) (-N/2+1:-1)]';
    c = c.*coef;
    c(abs(c)<1e-12) = 0; 
    z(idx(j,1):idx(j,2)) = ifft(c);
end