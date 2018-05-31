function plottimestep(z,idx,t,lambdas,extension,shear,j,curA,origA,fe,iters,modifs,U4,T4,S4,k4,rhok,uDomain,typeplot,comp_D,Def,zDomain,surfstruct)
%PLOTTIMESTEP
%Plots surfactant covered interface, domain and errors for timestep j.
%
%OBS. Ugly input. Will be fixed in a later version.
%
%:param: all information for plotting
%


if size(idx,1) == 1

    my_sfigure(1);

    rhomax = max(surfstruct.rho(idx(1,1):idx(1,2)));
    rhomin = min(surfstruct.rho(idx(1,1):idx(1,2)));

    %Plot the updated bubbles
    subplot(6,2,1:4);
    hold off
    cb = colorlinesplot(real(z),imag(z), ...
        surfstruct.rho,1,2,0.1,[],rhomin,rhomax,1);
    ylabel(cb,'Surfactant concentration')
    grid on
    axis equal
    hold on
    title(['$t=$ ',num2str(t), ' for $\lambda=$ ', num2str(lambdas), ' and $Q=$ ', num2str(extension), ', $G=$ ', num2str(shear)],'interpreter','latex');

    %The area error
    subplot(6,2,5);
    semilogy(j,abs(curA-origA)/origA,'*')
    hold on
    grid on
    title('Area error');

    %RK and GMRES steps
    subplot(6,2,6);
    title('GMRES (o), RK steps (*)');
    plot(j,fe/3,'*',j,iters/fe,'o',j,modifs/fe,'^')
    hold on
    grid on

    subplot(6,2,7);
    title('Npoints');
    plot(j,length(z),'*')
    grid on
    hold on

    %The number of quadrature modifications
    subplot(6,2,8);
    title('SQ modifications')
    semilogy(j,modifs/fe,'*')
    hold on
    grid on


    subplot(6,2,9)
    title('Conservation error surfactants, local(*) global (o)')
    semilogy(j,surfstruct.rhocons(end),'*')
    semilogy(j,surfstruct.rhoconsTOT(end),'o')
    grid on
    hold on

    subplot(6,2,10)
    hold off
    title('Surfactant concentration')
    hold on
    N = length(surfstruct.rho);
    alpha = linspace(0,2*pi,N+1)'; alpha = alpha(1:end-1);
    plot(alpha,surfstruct.rho,'-')
    xlabel('alpha')
    grid on


    subplot(6,2,11)
    N = length(surfstruct.rhok);
    kvec = (-N/2:N/2-1)';
    semilogy(kvec,abs(surfstruct.rhok)+eps,'.-')
    grid on
    title('Surfactants Fourier space')

    subplot(6,2,12)
    title('z fourier space')
    semilogy(kvec,abs(fftshift(1/N*fft(z))))
    hold on
    grid on

    if ~isempty(zDomain)
        sfigure(2)
        plot_domain(z,uDomain,typeplot,comp_D,extension)
    end

else
    my_sfigure(1);

    %Plot the updated bubbles
    subplot(5,2,1:4);
    rhomax = 0; rhomin = inf;
    for ji=1:size(idx,1)
        if max(surfstruct.rho(idx(ji,1):idx(ji,2))) > rhomax
            rhomax = max(surfstruct.rho(idx(ji,1):idx(ji,2)));
        end
        if min(surfstruct.rho(idx(ji,1):idx(ji,2))) < rhomin
            rhomin = min(surfstruct.rho(idx(ji,1):idx(ji,2)));
        end
    end
    hold off
    for ji=1:size(idx,1)
        cb = colorlinesplot(real(z(idx(ji,1):idx(ji,2))),imag(z(idx(ji,1):idx(ji,2))), ...
            surfstruct.rho(idx(ji,1):idx(ji,2)),1,[],0.1,[],rhomin,rhomax,1);
        grid on
        axis equal
        hold on
    end
    ylabel(cb,'Surfactant concentration');
    title(strcat('$t=$ ',num2str(t)));

    %The area error
    subplot(5,2,5);
    semilogy(j,abs(curA-origA)/origA,'*')
    hold on
    grid on
    title('Area error');

    %RK and GMRES steps
    subplot(5,2,6);
    plot(j,fe/3,'*',j,iters/fe,'o')
    hold on
    grid on
    title('GMRES (o), RK steps (*)');

    %Special quadrature modifications
    subplot(5,2,7);
    semilogy(j,modifs/fe,'*')
    hold on
    grid on
    title('SQ modifications');


    subplot(5,2,8);
    plot(j,length(z),'*')
    grid on
    hold on
    title('Npoints');

    if abs(rhomax) > 1e-12
        subplot(5,2,9);
        hold off
        for ji=1:size(idx,1)
            rho = [surfstruct.rho(idx(ji,1):idx(ji,2)); surfstruct.rho(idx(ji,1))];
            alpha = linspace(0,2*pi,length(rho))';
            plot(alpha,rho)
            grid on
            hold on
        end
        title('Surfactant concentration')
        ylim([0 rhomax])
        xlim([0 2*pi])

        subplot(5,2,10);
        for ji=1:size(idx,1)
            semilogy(j,surfstruct.rhocons(ji),'*')
            semilogy(j,surfstruct.rhoconsTOT(ji),'o')
            grid on
            hold on
        end
        title('Conservation error surfactants, local(*) global (o)')
    end


end

drawnow
