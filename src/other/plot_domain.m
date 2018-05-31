function plot_domain(z,uDomain,typeplot,comp_D,extension)
figure(99)
hold off

[zDomain,zplotOuter,zplotInner,Nouter] = comp_D(z);

switch typeplot
    case 'moving'
%         UDouter = reshape(uDomain(1:Nouter),size(zplotOuter));
%         pcolor(real(zplotOuter),imag(zplotOuter),log10(abs(UDouter)));
%         hold on
%         UDinner = reshape(uDomain(Nouter+1:end),size(zplotInner));
%         pcolor(real(zplotInner),imag(zplotInner),log10(abs(UDinner)));
%         axis([-3 3 -3 3])
%         
        UDouter = reshape(uDomain(1:Nouter),size(zplotOuter));
        pcolor(real(zplotOuter),imag(zplotOuter),abs(UDouter));
        hold on
        UDinner = reshape(uDomain(Nouter+1:end),size(zplotInner));
        pcolor(real(zplotInner),imag(zplotInner),abs(UDinner));
        axis([-3 3 -3 3])

%         figure(98)
%         hold off
%         ufarfield = extension*conj(zDomain(1:endOuter));
%         UFF = reshape(ufarfield,size(zplotOuter));
%         pcolor(real(zplotOuter),imag(zplotOuter),log10(abs(UDouter-UFF)));
%         shading flat
%         axis equal
%         hold on
%         plot(z,'k-')
%         colorbar
%         title('Difference to far field u')
        
    case 'box'
        UD = reshape(uDomain,size(zplot0));
        pcolor(real(zplot),imag(zplot),log10(abs(UD)))
        hold on
end
figure(99)
% quiver(real(zDomain),imag(zDomain),real(uDomain),imag(uDomain),'k','Autoscale','Off')
colorbar
shading flat
axis equal
plot(z,'k-')
title('u in domain')




end