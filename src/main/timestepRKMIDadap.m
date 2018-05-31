function [z,h,t,k2,U3,T2,S2,fe,iters,modifs,surf,dtinfo] = timestepRKMIDadap(z,t,dt,h,idx,W,bubble,omegaguess,beta,gamma,surf,vel_f,timetol,dtmax)
%TIMESTEPRKMIDADAP
%Goes from time t to t+dt with an adaptive Implicit-Explicit (IMEX)
%Runge Kutta method of order two. Adaptive both in space (z) and
%surfactants (rho)
%
% [z,h,t,k2,U3,T2,S2,fe,iters,modifs,surf] = timestepRKMIDadap(z,t,dt,h,idx,W,bubble,omegaguess,beta,gamma,surf,vel_f,timetol,dtmax)
%
%Returns:
% **z** -- new interface positions at time ts(j)
%
% **h** -- new timestep size, adapted after errors
%
% **k2,U3,T2,S2** -- velocities at new time t+dt
%
% **fe** -- number of Stokes evaluation from ts(j-1) to ts(j)
%
% **iters** -- number of GMRES iterations needed from ts(j-1) to ts(j)
%
% **modifs** -- number of special quadrature modifications from ts(j-1) to ts(j)
%
% **surf** -- struct containing all surfactant information at new time ts(j+1)
%
% **dtinfo** -- information regarding time step
%
%:param z: equidistant interface discretization
%:param t: current time
%:param dt: step size until return of function
%:param h: initial size of time step
%:param idx: index vector for multiple droplets
%:param W: quadrature weights
%:param bubble: bubble vector
%:param beta,gamma: viscosity parameters
%:param omegaguess: initial guess omega
%:param surf: struct containing surfactant information at time t=ts(j-1)
%:param vel_f: function to call for computing velocity
%:param timetol: adaptive timestepping tolerance
%:param dtmax: maximum time-step size
%

tend = t+dt;

fe = 0;
iters = 0;
modifs = 0;

dtinfo.nbrdt = 0;
dtinfo.nbrsuccess = 0;


curh = h;
while 1
    if h>dtmax
        h = dtmax;
%         disp(['Maxiumum size of timestep is: ' num2str(dtmax)]);
    end

    if t+h >= tend
        curh = h;
        h = tend-t;
    end

    %Evaluate things at time tn <- This is always done in the same way, for
    %all adaptive steps
    [k1,i1,m1,omegaguess,U1,T1,S1] = ...
        vel_f(z,idx,bubble,omegaguess,beta,gamma,surf.sigma);
    fe = fe + 1;
    iters = iters + i1;
    modifs = modifs + m1;
    frhok1 = comp_surf(surf.rhok,surf.salpha,surf.kappa,U1,T1,S1,idx); frhok1(1) = 0;

    while 1

        z1 = z + 0.5*h*k1;
        [zp1,zpp1] = fft_diff(z1,idx);

        salpha1 = zeros(size(surf.salpha)); kappa1 = zeros(size(z)); rhokhalf = zeros(size(surf.rhok)); rhohalf = zeros(size(surf.rho));
        for jidx=1:size(idx,1)
            I = idx(jidx,1):idx(jidx,2);
            N = length(z(I));
            kvec = (-N/2:N/2-1)';

            salpha1(jidx) = sum(W(I).*abs(zp1(I)))/2/pi;
            kappa1(I) = imag(zpp1(I)./zp1(I))./abs(zp1(I));

            %Evaluate at tn+1/2
            rhokhalf(I) = (surf.rhok(I) + 0.5*h*frhok1(I))./(1+0.5*h*kvec.^2/(surf.Pe(jidx)*salpha1(jidx)^2));
            rhokhalf(I(end)) = 0;
            rhohalf(I) = real(ifft(N*ifftshift(rhokhalf(I))));
        end

        sigmahalf = compute_surftension(rhohalf,surf.E,idx,surf.eqstate);
        [k2,i2,m2,omegaguess,U2,T2,S2] = ...
            vel_f(z1,idx,bubble,omegaguess,beta,gamma,sigmahalf);
        fe = fe +1;
        iters = iters + i2;
        modifs = modifs + m2;

        z2 = z + h*k2;
        [zp2,zpp2] = fft_diff(z2,idx);


        %Take one step with explicit Euler as first order method to compare
        %(position)
        zw = z + h*k1;
        relerrz = norm(z2 - zw)/norm(z2);


        frhok2 = comp_surf(rhokhalf,salpha1,kappa1,U2,T2,S2,idx); frhok2(1)=0;

        rhok2 = zeros(size(surf.rhok)); rho2 = zeros(size(surf.rho));
        salpha2 = zeros(size(salpha1)); kappa2 = zeros(size(kappa1));
        rhoks = zeros(size(surf.rhok));
        for jidx=1:size(idx,1)
            I = idx(jidx,1):idx(jidx,2);
            N = length(z(I));

            salpha2(jidx) = sum(W(I).*abs(zp2(I)))/2/pi;
            kappa2(I) = imag(zpp2(I)./zp2(I))./abs(zp2(I));

            kvec = (-N/2:N/2-1)';
            rhok2(I) = surf.rhok(I) + h*(-kvec.^2).*rhokhalf(I)/(surf.Pe(jidx)*salpha1(jidx).^2) + h*frhok2(I);
            rhok2(I(end)) = 0;
            rho2(I) = real(ifft(N*ifftshift(rhok2(I))));



            rhoks(I) = (surf.rhok(I) + h*frhok1(I))./(1+h*kvec.^2/(surf.Pe(jidx)*salpha2(jidx)^2));
            rhoks(I(end)) = 0;
        end
        sigma2 = compute_surftension(rho2,surf.E,idx,surf.eqstate);

        % This velocity is only needed as for adaptive-rho comp. with IMEX1 below
        [~,~,~,omegaguess,U3,T3,S3] = ...
            vel_f(z2,idx,bubble,omegaguess,beta,gamma,sigma2);
%        fe = fe + 1; iters = iters + i2; modifs = modifs + m2;
%         frhok3 = comp_surf(rhoks,salpha2,kappa2,U3,T3,S3,idx);
%         rhokw = zeros(size(surf.rhok));
%         for jidx=1:size(idx,1)
%             I = idx(jidx,1):idx(jidx,2);
%             N = length(z(I));
%             kvec = (-N/2:N/2-1)';
%             rhokw(I) = surf.rhok(I) + h*(-kvec.^2).*rhoks(I)/(surf.Pe(jidx)*salpha2(jidx)^2) + h*frhok3(I);
%             rhokw(I(end)) = 0;
%         end

        if surf.rhomass0(1) ~= 0
            %Adaptivity in rho with local conservation
            relerrrho = 0;
            for jidx=1:size(idx,1)
                I = idx(jidx,1):idx(jidx,2);
                mass2 = abs(sum(rho2(I).*salpha2(jidx).*W(I)));
                relerrrho = max(relerrrho,abs(mass2-surf.rhomass(jidx,end))/abs(surf.rhomass(jidx,end)));
            end

%             %Adaptivity in rho with IMEX1
%            relerrrho = norm(rhok2-rhokw)/norm(rhok2);

        else
            relerrrho = 0;
        end
        relerr = max(relerrrho,relerrz);
%         relerr = relerrrho;

        dtinfo.nbrdt = dtinfo.nbrdt + 1;
        if timetol ~= 0
            if relerr > timetol
                disp(['h = ' num2str(h) '. Fail, decrease time step'])
                h = h*(0.9*timetol/relerr)^(1/2);
            else
                disp(['h = ' num2str(h) '. Success'])
                dtinfo.nbrsuccess = dtinfo.nbrsuccess + 1;
                break
            end
        else
            disp(['No adaptive timestep, h = ' num2str(h)])
            break;
        end

    end
    z = z2;
    [zp,zpp] = fft_diff(z,idx);

    % FILTER
    z = fft_smooth(z,idx,0.1);
    rho2 = real(fft_smooth(rho2,idx,0.1));
    for k=1:size(idx,1)
        I = idx(k,1):idx(k,2);
        N = length(I);
        rhok2(I) = 1/N*fftshift(fft(rho2(I)));
    end
    rhok2(abs(rhok2) < 1e-12) = 0; %KRASNY

    for jidx=1:size(idx,1)
        I = idx(jidx,1):idx(jidx,2);
        surf.salpha(jidx) = sum(W(I).*abs(zp(I)))/2/pi;
        surf.kappa(I) = imag(zpp(I)./zp(I))./abs(zp(I));
        surf.rhok(I) = rhok2(I);
        surf.rho(I) = rho2(I);

        surf.rhomass(jidx,2) = abs(sum(surf.rho(I).*surf.salpha(jidx).*W(I)));
        e = abs(surf.rhomass(jidx,2)-surf.rhomass(jidx,1))/surf.rhomass(jidx,1);
        tmp = max(surf.rhocons(jidx),e);
        surf.rhocons(jidx) = tmp;
        surf.rhoconsTOT(jidx) = abs(surf.rhomass(jidx,2)-surf.rhomass0(jidx))/surf.rhomass0(jidx);
    end
    surf.sigma = compute_surftension(surf.rho,surf.E,idx,surf.eqstate);

    t = t+h;
    if abs(t-tend) < 1e-13
        break;
    end

    if timetol~=0
        h = h*(0.9*timetol/relerr)^(1/2);
    end
end
h = curh;
end
