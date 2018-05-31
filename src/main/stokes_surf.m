function [final] = stokes_surf(ipF,filename,ifPlot)
%STOKES - SURF
%GL-Trapezoidal hybrid for solving 2D-Stokes eq's on multiple interfaces
%having different lambdas.
%
%Version 2 - added insoluble surfactants.
%
%Date of original version (2.0) : 2015-11-16
%
%Updated: 2016-12-01.
%
% [final] = stokes_surf(ipF,filename,ifPlot)
%
%Returns:
%  **final** -- struct containing simulation data at final timestep
%
%:param ipF: inputstructfile from which to read set parameters
%:param filename: save path for output data
%:param ifPlot: boolean if plotting during simulation or not
%

if nargin < 3
    ifPlot = 0;
    if nargin < 2
        filename = 'output/tmp';
        if nargin < 1
            error('ERROR: no input file. Simulation stopped.');
        end
    end
end


[globaldata,timedisc,drops,surfact] = readInput(ipF);

D = D16init;

warning on

disp(['Number of points: ' num2str(length(drops.z))]);
disp(['Viscosity ratio set to: ' num2str(drops.lambdas')]);
disp(['Capillary number set to: ' num2str(globaldata.extension)]);
if globaldata.breakwhensteady
    disp('Stop simulation at steady state.')
end
disp(['Run from time ' num2str(timedisc.tSpan(1)) ' to time ' num2str(timedisc.tSpan(2))])

if ifPlot
    scrsz = get(0,'ScreenSize'); %MOVIE
    my_sfigure(1,'Name','Bubble chamber','Position',[1 scrsz(4)/2 2*scrsz(3)/4 3*scrsz(4)/4]); %MOVIE

    plottimestep(drops.z,drops.idx,timedisc.t,drops.lambdas,globaldata.extension,globaldata.shear,0,drops.area,globaldata.origA,0,0,0,0,0,0,0,[],[],[],[],0,[],surfact)

end

mastertime = tic;


%The function computeVelocity computes the resulting velocity of
%a boundary consisting of the points z.
vel_f = @(z,idx,bubble,omegaguess,beta,gamma,sigma)  ...
    computeVelocity(z,idx,bubble,D,beta,gamma,globaldata.shear, ...
    globaldata.extension,omegaguess,sigma);

% - - - - - - - - - - - TIME STEP - - - - - - - - - - - - - - - - - - - -
omegaguess = zeros(2*length(drops.z),1);
sim_stalled = 0;
steadyBreak = 0;
feTOT = 0; %total number of function evaluations
dtTOTsuccess = 0; %total number of time-steps taken
dtTOT = 0;

% Starting time
t = timedisc.t;

drops.oldarea = drops.area; %OBS WRONG FIRST GO

for j = timedisc.startj:timedisc.nT+1


    % - - - - - - - - - - SPECIFY TIMESTEP - - - - - - - - - - - - - - - -
    [z,~,t,k4,U4,T4,S4,fe,iters,modifs,surfact,dtinfo] = ...
        timestepRKMIDadap(drops.z,t,timedisc.dt,timedisc.h,drops.idx,drops.W,drops.bubble,omegaguess, ...
        drops.beta,drops.gamma,surfact,vel_f,globaldata.timeTol,globaldata.time.dtmax);
    drops.z = z; %Remember to also fix area!
    feTOT = feTOT + fe;
    dtTOT = dtTOT + dtinfo.nbrdt;
    dtTOTsuccess = dtTOTsuccess + dtinfo.nbrsuccess;


    % Check if we should filter:
    if floor((j-1)/timedisc.Nf)-ceil((j-1)/timedisc.Nf) == 0
        drops.z = fft_smooth(z,drops.idx,0.1);
        if globaldata.surfact.rhoinit ~= 0
            surfact.rho = real(fft_smooth(surfact.rho,drops.idx,0.1));
            for k=1:size(drops.idx,1)
                I = drops.idx(k,1):drops.idx(k,2);
                N = length(I);
                surfact.rhok(I) = 1/N*fftshift(fft(surfact.rho(I)));
            end
            surfact.rhok(abs(surfact.rhok) < 1e-12) = 0; %KRASNY
            surfact.sigma = compute_surftension(surfact.rho,surfact.E,drops.idx,surfact.eqstate);
        end
    end

    %If we want to break when the bubbles are approximately circular, we
    %need to check if this is indeed the case.
    if globaldata.breakwhencircle ~= 0
        maxfactor = -1e78;
        for k = 1:size(drops.idx,1)
            curz = drops.z(drops.idx(k,1):drops.idx(k,2));
            deviation = abs(curz-mean(curz));
            factor = abs(1-max(deviation)/mean(deviation));
            if factor > maxfactor
                maxfactor =factor;
            end
        end
        if maxfactor < globaldata.breakwhencircle
            disp('Drops are circular enough')
            break;
            %         else %RIKARDS THING. OBSOLETE?
            %             timedisc.ts(j+1) = timedisc.ts(j)+suph;
        end
    end

    % If we want to break for steady state we need to check when it's
    % reached.
    if globaldata.breakwhensteady
        Unorm = norm(U4,inf);
        if Unorm < globaldata.breakwhensteady
            disp(['Steady state reached! unorm = ' num2str(Unorm)])
            steadyBreak = 1;
            if ~isempty(zDomain)
                sfigure(2)
                plot_domain(drops.z,uDomain,typeplot,comp_D,globaldata.extension)
            end
            break;
        end
    end

    % Print out area and conservation error
    zp = fft_diff(drops.z,drops.idx);
    curA = abs(sum(drops.W.*real(drops.z).*imag(zp)));
    curL = zeros(size(drops.idx,1),1);
    %The maximum arclength ratio and the bubble boundary lengths
    curmax = -1;
    Def = zeros(drops.nbrDrops,1);
    for k = 1:size(drops.idx,1)
        I = drops.idx(k,1):drops.idx(k,2);
        curL(k) = sum(drops.W(I).*abs(zp(I)));
        az = abs(diff(drops.z(I)));
        if max(az)/min(az) > curmax
            curmax = max(az)/min(az);
        end
        equalarc = curL(k)/(2*pi);
        % Calculate distances to obtain D
        Def(k) = compDeform(drops.z(drops.idx(k,1):drops.idx(k,2)));
    end
    drops.area = curA;
    % Switch surfactant mass over for previous time-step
    surfact.rhomass(:,1) = surfact.rhomass(:,2);

    %We are now done with the current step. Check the npanel/length ratios
    %for each panel. If the boundary of the bubble is too long, shrink the
    %number of discretization points describing it. If too short, extend.
    if globaldata.doadaptivity == 1
        if floor((j-1)/timedisc.Na)-ceil((j-1)/timedisc.Na) == 0
        disp('Check spatial adaptivity..');

            nbrDrops = size(drops.idx,1);

            for k = 1:nbrDrops

                idx = drops.idx;
                oldidx = idx;

                newN = ceil(globaldata.origq(k)*curL(k))+1;
                if (idx(k,2)-idx(k,1)+1)/16 > newN
                    sdiff = (idx(k,2)-idx(k,1)+1) - newN*16;
                    newz = shrink(drops.z(idx(k,1):idx(k,2)),newN*16);
                    drops.z = [drops.z(1:idx(k,1)-1);newz;drops.z(idx(k,2)+1:end)];

                    %Extend rho to fit new discretization
                    newrho = real(shrink(surfact.rho(idx(k,1):idx(k,2)),newN*16));
                    surfact.rho = [surfact.rho(1:idx(k,1)-1); newrho; surfact.rho(idx(k,2)+1:end)];

                    drops.bubble(2*idx(k,1)-1 + (0:2*sdiff-1)) = [];
                    drops.lambda(2*idx(k,1)-1 + (0:2*sdiff-1)) = [];
                    drops.beta(2*idx(k,1)-1 + (0:2*sdiff-1)) = [];
                    drops.gamma(2*idx(k,1)-1 + (0:2*sdiff-1)) = [];

                    drops.W = [drops.W(1:idx(k,1)-1);2*pi/newN/16*ones(length(newz),1);drops.W(idx(k,2)+1:end)];

                    idx(k,2) = idx(k,2)-sdiff;
                    for l = k+1:size(idx,1)
                        idx(l,:) =idx(l,:) - sdiff;
                    end

                    %Update kappa in surf
                    [zp,zpp] = fft_diff(drops.z(idx(k,1):idx(k,2)),[1 length(drops.z(idx(k,1):idx(k,2)))]);
                    newkappa = imag(zpp./zp)./abs(zp);
                    surfact.kappa = [surfact.kappa(1:oldidx(k,1)-1); newkappa; surfact.kappa(oldidx(k,2)+1:end)];

                    drops.idx = idx;

                    drops.Nps(k) = newN;

                    disp('Number of discretization points decreased');
                elseif (idx(k,2)-idx(k,1)+1)/16 < newN
                    sdiff = newN*16-(idx(k,2)-idx(k,1)+1);
                    newz = shrink(drops.z(idx(k,1):idx(k,2)),newN*16);
                    drops.z = [drops.z(1:idx(k,1)-1);newz;drops.z(idx(k,2)+1:end)];

                    %Extend rho to fit new discretization
                    newrho = real(shrink(surfact.rho(idx(k,1):idx(k,2)),newN*16));
                    surfact.rho = [surfact.rho(1:idx(k,1)-1); newrho; surfact.rho(idx(k,2)+1:end)];

                    drops.lambda = [drops.lambda(1:2*idx(k,1)-2);drops.lambda(2*idx(k,1)-1)*ones(2*sdiff,1);drops.lambda(2*idx(k,1)-1:end)];
                    drops.bubble = [drops.bubble(1:2*idx(k,1)-2);drops.bubble(2*idx(k,1)-1)*ones(2*sdiff,1);drops.bubble(2*idx(k,1)-1:end)];
                    drops.beta = [drops.beta(1:2*idx(k,1)-2);drops.beta(2*idx(k,1)-1)*ones(2*sdiff,1);drops.beta(2*idx(k,1)-1:end)];
                    drops.gamma = [drops.gamma(1:2*idx(k,1)-2);drops.gamma(2*idx(k,1)-1)*ones(2*sdiff,1);drops.gamma(2*idx(k,1)-1:end)];

                    drops.W = [drops.W(1:idx(k,1)-1);2*pi/newN/16*ones(length(newz),1);drops.W(idx(k,2)+1:end)];
                    idx(k,2) = idx(k,2)+sdiff;
                    for l = k+1:size(idx,1)
                        idx(l,:) =idx(l,:) + sdiff;
                    end
                    %Update kappa in surf
                    [zp,zpp] = fft_diff(drops.z(idx(k,1):idx(k,2)),[1 length(drops.z(idx(k,1):idx(k,2)))]);
                    newkappa = imag(zpp./zp)./abs(zp);
                    surfact.kappa = [surfact.kappa(1:oldidx(k,1)-1); newkappa; surfact.kappa(oldidx(k,2)+1:end)];

                    drops.idx = idx;

                    drops.Nps(k) = newN;

                    disp('Number of discretization points increased');
                end
            end
        end
        drops.z = fft_smooth(drops.z,drops.idx,0.1);
        surfact.rho = real(fft_smooth(surfact.rho,drops.idx,0.1));
        for k=1:size(drops.idx,1)
            I = drops.idx(k,1):drops.idx(k,2);
            N = length(I);
            surfact.rhok(I) = 1/N*fftshift(fft(surfact.rho(I)));
        end
        surfact.rhok(abs(surfact.rhok) < 1e-12) = 0; %KRASNY
        surfact.sigma = compute_surftension(surfact.rho,surfact.E,drops.idx,surfact.eqstate);
    end

    % Plot data
    if ifPlot && floor((j-1)/timedisc.Np)-ceil((j-1)/timedisc.Np)==0
        disp('Plot!')
        plottimestep(drops.z,drops.idx,t,drops.lambdas,globaldata.extension,...
            globaldata.shear,j,curA,globaldata.origA,fe,iters,modifs,U4,T4,S4,k4, ...
            [],[],[],[],Def,[],surfact);
    end

    % Save data
    if floor((j-1)/timedisc.Ns)-ceil((j-1)/timedisc.Ns) == 0
        disp('Save data...');
        simdata.drops = drops;
        simdata.surfact = surfact;
        simdata.t = t;
        simdata.globaldata = globaldata;
        simdata.nbrsuccess_dt = dtTOTsuccess;
        simdata.nbrtot_dt = dtTOT;
        filename_tmp = [filename '_t' num2str(t)];
        filename_tmp(filename_tmp == '.') = 'p';
        save(filename_tmp,'simdata');
        disp(['Data for time ' num2str(t) ' saved in file: ' filename_tmp])
    end

    if floor((j-1)/timedisc.Np)-ceil((j-1)/timedisc.Np)==0
        disp('Print info...')
        disp('')
        disp(['time = ',num2str(t)])
        if globaldata.surfact.rhoinit ~= 0
            disp(['Conservation error surf: ' num2str(max(surfact.rhocons)) ' (local)'])
            disp(['Conservation error surf: ' num2str(max(surfact.rhoconsTOT)) ' (global)'])
        end
        disp(['Area error: ' num2str(abs(drops.area-drops.oldarea)/globaldata.origA) ' (local)'])
        disp(['Area error: ' num2str(abs(drops.area-globaldata.origA)/globaldata.origA) ' (global)'])
        disp(['Deformation: ' num2str(max(Def))])
        if globaldata.breakwhensteady > 0
            disp(['Normal velocity: ' num2str(Unorm)])
        end
        disp('');

    end
        drops.oldarea = drops.area;
end

Ae = abs(curA-globaldata.origA)/globaldata.origA;
err.Ae = Ae;
err.rhoconsTot = surfact.rhoconsTOT;
err.rhocons = surfact.rhocons;

totaltime = toc(mastertime);

disp('Done!');
disp(['Area error: ' num2str(Ae)]);
disp(['Conservation error: ' num2str(max(err.rhoconsTot))]);
disp(['Total time: ' num2str(totaltime)]);

% Save data at the end!
simdata.drops = drops;
simdata.surfact = surfact;
simdata.t = t;
simdata.globaldata = globaldata;
simdata.nbrsuccess_dt = dtTOTsuccess;
simdata.nbrtot_dt = dtTOT;
filename_tmp = [filename '_t' num2str(t)];
filename_tmp(filename_tmp == '.') = 'p';
save(filename_tmp,'simdata');
disp(['Data for time ' num2str(t) ' saved in file: ' filename_tmp])

final.drops = drops;
final.D = Def;
final.t = t;
final.time = totaltime;
final.surfact = surfact;
final.err = err;
final.stokes = feTOT;
final.steadybreak = steadyBreak;
final.stalled = sim_stalled;
final.nbrsuccess_dt = dtTOTsuccess;
final.nbrtot_dt = dtTOT;
