function [zDomain,zplotOuter,zplotInner,NOuter] = comp_domain(z,typeDisc)
 
NOuter = 0;
switch typeDisc
    case 'circle'
        indjump = 10;
        R = [linspace(0,0.99,30)'; linspace(0.01,2,40)'];
        ztest = z(1:indjump:end);
        T = angle(ztest); T = [T; T(1)];
        [Rplot,Tplot] = meshgrid(R,T);
        zDomain = Rplot(:).*exp(1i*Tplot(:));
        zplotOuter = Rplot.*exp(1i*Tplot);
        zplotInner = [];
        
    case 'box'
        % Box in which to put domain points
        % NB some points that we distribute might come too close to the boundary,
        % as the boundary moves. so we might get large errors very close to the
        % interface... Working on that
        x0 = -5; x1 = 5;
        y0 = -5; y1 = 5;
        X = linspace(x0,x1,100)';
        Y = linspace(y0,y1,100)';
        % X = 0.1; Y = 0.1;
        [xx, yy] = meshgrid(X,Y);
        zplotOuter = xx + 1i*yy;
        zDomain = reshape(zplot,length(zplot)^2,1);
        Rplot = zeros(length(zplot));
        zplotInner = [];
        
    case 'moving'
        indjump = 5;
        ztest = z(1:indjump:end);
        %Outer domain
        Rmax = 5; 
        Nr = 500;
        alpha = linspace(0,2*pi,length(ztest)+1)'; %alpha = alpha(1:end-1);
        zouter = (cos(alpha) + 1i*sin(alpha))*Rmax;
        [X,Y] = GenericPolarGrid([[real(ztest); real(ztest(1))] [imag(ztest);imag(ztest(1))]],[real(zouter) imag(zouter)],Nr);
        zplotOuter = X + 1i*Y;
        Rplot = X;
        zDomain1 = zplotOuter(:);
        %Inner domain
        Nr = 50;
        zorigo = zeros(size(ztest));
        [X2,Y2] = GenericPolarGrid([[real(ztest); real(ztest(1))] [imag(ztest); imag(ztest(1))]],[[real(zorigo); real(zorigo(1))] [imag(zorigo);imag(zorigo(1))]],Nr);
        zplotInner = X2 + 1i*Y2;
        zDomain2 = zplotInner(:);
        
        zDomain = zDomain1;
%         zDomain = [zDomain1; zDomain2];
        NOuter = length(zDomain1);
end
end