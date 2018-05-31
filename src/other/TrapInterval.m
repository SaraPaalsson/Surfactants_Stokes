function [T,W] = TrapInterval(a,b,N)
T = linspace(a,b,N+1)';
T = T(1:end-1);
W = ones(N,1)*(b-a)/N;