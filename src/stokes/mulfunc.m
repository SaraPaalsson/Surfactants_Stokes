function y = mulfunc(x,z,zp,z_sc,zp_sc,zpp,W,beta,zc,cidx,imodifs,didx) 
%MULFUNC
%Computes matrix-vector multiplication (through FMM) for solving Ax=b for
%x. Called by gmres.  
%
%Returns:
%  **y** -- computed matrix-vector product
%
%:param x: vector for multiplying (i.e. old iteration of complex density
%           gmres)
%:param  z,zp,zpp: discretization points and their derivatives
%:param z_sc,zp_sc: scaled discretization points and derivatives
%:param W: Gauss-Legendre quadrature weights
%:param beta: parameter representing viscosity ratio of problem
%:param zc: centre points of droplets
%:param cidx: index table for using special quadrature
%:param imodifs: how to modify solution with special quadrature
%:param didx: index table
%



[~,y] = mex_M1M4(z_sc,zp_sc,W.*x/pi,1e-13,8,4);

y = -1i*beta.*y + beta.*W.*x.*imag(zpp./zp)/2/pi + beta.*W.*conj(x).*imag(zpp.*conj(zp))./(conj(zp).^2)/2/pi ;

imodifs = -imodifs;
mex_applymodifs(y,x,imodifs,cidx);

%For drops, viscous outside
zc = zeros(size(didx,1),1); pk = zeros(size(didx,1),1);
for idrops = 1:size(didx,1)
    zc(idrops) = sum(z(didx(idrops,1):didx(idrops,2)))/length(z(didx(idrops,1):didx(idrops,2)));
    pk(idrops) = 2/pi*real(sum(W(didx(idrops,1):didx(idrops,2)).* ...
        x(didx(idrops,1):didx(idrops,2)).*zp(didx(idrops,1):didx(idrops,2))./ ...
        ((z(didx(idrops,1):didx(idrops,2))-zc(idrops)).^2)));
end

pinf = -pk(1);
pk = pinf + pk;

for idrops = 1:size(didx,1)
    
    if idrops == 1
        j1 = 1; j2 = didx(idrops,2);
    else
        j1 = didx(idrops-1,2)+1;
        j2 = didx(idrops,2);
    end
    Ht = sum(W(didx(idrops,1):didx(idrops,2)).*...
        abs(zp(didx(idrops,1):didx(idrops,2))).*x(didx(idrops,1):didx(idrops,2)));
    
    y(j1:j2) = y(j1:j2) - 0.5*beta(j1:j2).*1i*(pinf-pk(idrops)).*z(j1:j2) ...
        + beta(j1:j2).*Ht;
end




