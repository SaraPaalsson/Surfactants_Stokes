function c = bandstop(N,proc,a)
f1 = floor(N/2-proc*N);
f2 = ceil(N/2+proc*N);
t = (1:N)';
c1 = 1- 1./(1+exp(-a*(t-f1)));
c2 =1-1./(1+exp(-a*(N-t-f1)));
c = c1+c2;

end