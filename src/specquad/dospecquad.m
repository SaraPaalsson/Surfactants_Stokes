function [velocity,nmodifs] = dospecquad(velocity,omega,z,zp,W,pe,bubble,idx)
%DOSPECQUAD
%Matlab version of the function that computes special quadrature
%modifications. 
%
%OBS: Slow.
%
%  [velocity,nmodifs] = dospecquad(velocity,omega,z,zp,W,pe,bubble,idx)
%
%Returns:
%  **velocity** -- updated velocity
%
%  **nmodifs** -- number of modifications made by special quadrature
%
%:param z,zp: interface discretization points and derivatives
%:param W: quadrature weights
%:param idx: index vector for multiple droplets
%:param velocity: computed velocity to be modified
%:param pe: panel edges
%:param bubble: bubble vector
%:param omega: complex density
%

p = zeros(16,1);
q = zeros(16,1);

xmax = max(real(z));
xmin = min(real(z));
ymax = max(imag(z));
ymin = min(imag(z));

NPanels = length(z)/16;

pe2 = [];
ptr = 0;
for j = 1:size(idx,1)
    np = (idx(j,2)-idx(j,1)+1)/8;
    pe2 = [pe2;pe(j-1 + ptr + (1:np)) pe(j-1 + ptr + (2:np+1))];
    ptr = ptr + np;
end

panel_lengths = abs(pe2(:,2)-pe2(:,1));

meanlen = mean(panel_lengths);

XBoxes = ceil((xmax-xmin)/meanlen);
YBoxes = ceil((ymax-ymin)/meanlen);

grid = cell(XBoxes,YBoxes);

for n = 1:NPanels
    z1 = pe2(n,2);
    z2 = pe2(n,1);
    mid = (z1+z2)/2;
    zrel = (mid -(xmin+1i*ymin))/meanlen;
    midx = floor(real(zrel))+1;
    midy = floor(imag(zrel))+1;
    %Danger here
    radius = ceil(panel_lengths(n)/meanlen);
    for x = midx-radius:midx+radius
        for y = midy-radius:midy+radius
            if x > 0 && x <= XBoxes && y > 0 && y <= YBoxes
                grid{x,y} = [grid{x,y} n];
            end
        end
    end
end

for x=1:XBoxes
    for y=1:YBoxes
        grid{x,y} = unique(grid{x,y});
    end
end

nmodifs = 0;
%Special quadrature
for k = 1:length(z)
    zrel = (z(k) -(xmin+1i*ymin))/meanlen;
    midx = floor(real(zrel))+1;
    midy = floor(imag(zrel))+1;
    
    vk2 = 0;
    for j = grid{midx,midy}
        
        
        b1 = bubble((j-1)*16+1);
        mid = (pe(j+1+b1)+pe(j+b1))/2;
        len = pe(j+1+b1)-pe(j+b1);
        
        
        if abs(z(k)-mid) < abs(len) && ...
                (bubble(k) ~= bubble((j-1)*16+1) || (abs(j-floor((k-1)/16)-1) > 1 && abs(j-floor((k-1)/16)-1) ~= (idx(b1+1,2)-idx(b1+1,1)+1)/8-1))
            
            pnlidx = (j-1)*16+(1:16);
            nz = 2*(z(k)-mid)/len;
            nzpan = 2*(z(pnlidx)-mid)/len;
            p(1) = log(1-nz)-log(-1-nz);
            q(1) = -1/(1-nz)-1/(1+nz);
            old = W(pnlidx).*zp(pnlidx)./(z(pnlidx)-z(k));
            if imag(nz) > 0 && real(nz) > -1 && real(nz) < 1
                
                if any(imag(nzpan) > imag(nz))
                    %Expensive but accurate test.
                    cc = vandernewtonT(real(nzpan),imag(nzpan));
                    vv2 = cc(end);
                    for kk = 15:-1:1
                        vv2 = vv2*real(nz) + cc(kk);
                    end
                    if vv > imag(nz)
                        p(1) = p(1) - 2i*pi;
                    end
                end
                
            end
            
            
            if abs(p(1)-sum(old))>1e-14
                
                for kk = 2:16
                    p(kk) = nz*p(kk-1) + (1-(-1)^(kk-1))/(kk-1);
                    q(kk) = nz*q(kk-1) + p(kk-1);
                end
                
                pp = mex_vandernewton(nzpan,p);
                qq = mex_vandernewton(nzpan,q)*2/len;
                
                qold = old./(z(pnlidx)-z(k));
                
                M1mod = omega(k)*real(p(1)-sum(old));
                M1realmod = real(pp-old).'*omega(pnlidx);
                M5mod = conj((pp-old)./zp(pnlidx)).*zp(pnlidx) - conj(qq-qold).*(z(pnlidx)-z(k));
                vk2 = vk2 + (M1mod-M1realmod+M5mod.'*conj(omega(pnlidx))/2)/pi;
                nmodifs = nmodifs + 1;
            end
            
            
        end
    end
    velocity(k) = velocity(k) + vk2;
end