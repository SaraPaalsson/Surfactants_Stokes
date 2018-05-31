tol = 1e-10;
N = 10;


b = rand(N,1);
A = rand(N,N) + 1i*rand(N,N);

%% Test 1: compare to MATLABs gmres
xmatlab = gmres(@(x) gmres_KX(x,A,1),b);
xours = mygmresg_el(b,0,1e-10,1000,0,@(x,varargin) gmres_KX(x,A,2));

assert(max(abs(xmatlab-xours)) <= tol, 'GMRES ok')