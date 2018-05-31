function Dmin = findClosestDist(z,idx)
% Compute the closest distance between drops described by z array

N = size(idx,1); %Number drops

% Plot - take away later
figure(999); clf
for k=1:N
    plot(z(idx(k,1):idx(k,2)),'k')
    hold on
end
axis equal

% Go through all drops and find closest distance between drop k and all
% other drops
Dminvec = zeros(N,1);
zmin = zeros(N,2);
for k=1:N 
    if k>1
        zrest = z(1:idx(k-1,2));
    else
        zrest = [];
    end
    if k<N
        zrest = [zrest; z(idx(k+1,1):end)];
    end
    zk = z(idx(k,1):idx(k,2));
   
%     plot(zk,'r')
    
    X = [real(zrest) imag(zrest)];
    Y = [real(zk) imag(zk)];
    
    [Didxk,D] = knnsearch(X,Y);
    [Dmink,Imink] = min(D);
    
    zk_min = zk(Imink);
    zrest_min = zrest(Didxk(Imink));
    
    zmin(k,1) = zk_min;
    zmin(k,2) = zrest_min;
    
%     plot(zk_min,'*');
%     plot(zrest_min,'o');
    
    Dminvec(k) = Dmink;
   
%     pause(0.5)
%     plot(zk,'k')
end

[Dmin,Imin] = min(Dminvec);

figure(999);
plot(zmin(Imin,1),'*')
plot(zmin(Imin,2),'o')

end