function [globaldata,timedisc,drops,surfact] = ...
    readInput(ipF)
%readInput
%Loads inputstruct file ipF where all parameters for running stokes_surf are
%set. Computes discretization and sets all initial data for our
%simulation.
%
% [globaldata,timedisc,drops,surfact,datastruct] = readInput(ipF)
%
%Returns:
% a collection of structs containing all parameters, as well as initial discretization of space and surfactants
%
%:param ipF: inputstructfile from which to read set parameters
%


run (ipF)

if isfield(inputstruct,'isrestartfile')
    isrestart = inputstruct.isrestartfile;
    restartfile = inputstruct.restartfile;
else
    isrestart = 0;
end

if isrestart
    load(restartfile)
    drops = simdata.drops;
    surfact = simdata.surfact;
    globaldata = simdata.globaldata;
    globaldata.isrestart = isrestart;
    globaldata.restartfile = restartfile;
end

% Read in time discretization
timedisc.tSpan = inputstruct.timedisc.Tint; %time-interval

timedisc.nT = inputstruct.timedisc.Nt; %number of non-adaptive time-steps with adaptive substeps
timedisc.dtmax = inputstruct.timedisc.dtmax;

timedisc.hdivisor = inputstruct.timedisc.hdivisor;  %guess on size of substeps


% Check for parameters of save-intervals and adaptivity checks
if isfield(inputstruct,'actionintervals')
    timedisc.Ns = inputstruct.actionintervals.Ns;
    timedisc.Na = inputstruct.actionintervals.Na;
    timedisc.Nf = inputstruct.actionintervals.Nf; %Obosolete?
    timedisc.Np = inputstruct.actionintervals.Np;
else
    timedisc.Ns = timedisc.nT;
    timedisc.Na = timedisc.nT;
    timedisc.Nf = timedisc.nT;
    timedisc.Np = timedisc.nT;
end

% Compute time vector and timestep h
timedisc.startj = 2;
ts = linspace(timedisc.tSpan(1),timedisc.tSpan(2),timedisc.nT+1)';
timedisc.t = ts(1);
timedisc.h = (ts(2)-ts(1))/timedisc.hdivisor;
timedisc.dt = inputstruct.timedisc.dt;

% Save information for global struct
timediscglobal.tSpan = timedisc.tSpan;
timediscglobal.nT = timedisc.nT;
timediscglobal.dtmax = timedisc.dtmax;
timediscglobal.Ns = timedisc.Ns;
timediscglobal.Na = timedisc.Na;
timediscglobal.Nf = timedisc.Nf;

globaldata.time = timediscglobal;
globaldata.breakwhensteady = inputstruct.steadyTol;

% Read in globaldata
globaldata.timeTol = inputstruct.timeTol;
globaldata.breakwhencircle = inputstruct.circleTol;
globaldata.breakwhensteady = inputstruct.steadyTol;
globaldata.doadaptivity = inputstruct.doadaptivity;
globaldata.extension = inputstruct.farfield.Q;
globaldata.shear = inputstruct.farfield.G;


if ~isrestart

    % Read in globaldata
    globaldata.timeTol = inputstruct.timeTol;
    globaldata.breakwhencircle = inputstruct.circleTol;
    globaldata.breakwhensteady = inputstruct.steadyTol;
    globaldata.doadaptivity = inputstruct.doadaptivity;
    globaldata.extension = inputstruct.farfield.Q;
    globaldata.shear = inputstruct.farfield.G;
    globaldata.time = timediscglobal;

    globaldata.surfact.eqstate = inputstruct.surfact{1}.interf.eqstate;


    % Read in drop-info
    drops.nbrDrops = inputstruct.nbrDrops;
    drops.pointsPerPanel = inputstruct.drops{1}.pointsperpanel; %always same

    zFuncs_mine = cell(drops.nbrDrops,1);

    % Compute initial interface discretization
    for i = 1:drops.nbrDrops
        switch inputstruct.drops{i}.shape
            case 'ellipse'
                aspect = inputstruct.drops{i}.aspect;
                a = sqrt(aspect);
                b = 1/sqrt(aspect);
                c = inputstruct.drops{i}.centre;

                if isfield(inputstruct.drops{i},'theta')
                    theta = inputstruct.drops{i}.theta;
                else
                    theta = 0;
                end

                zFuncs_mine{i} = @(T) ellipse_func(a,b,c,theta,T);
            case 'flower'
                aspect = inputstruct.drops{i}.aspect;
                a = sqrt(aspect);
                b = 1/sqrt(aspect);
                c = inputstruct.drops{i}.centre;

                try
                    thet = inputstruct.drops{i}.theta;
                catch
                    thet = 0;
                end
                zFuncs_mine{i} = @(T) flower_func(a,inputstruct.drops{i}.scale,thet,inputstruct.drops{i}.Narms,c,T);

            case 'swissroll'
                load 'spiralnew'
                load 'swissroll'
                zFuncs_mine{i} = @(T) mex_evalspline(c,T);

            case 'swissroll_ellipse'
                load 'spiralnew'
                load 'swissroll'
%                 zFuncs_mine{i} = zFuncs{i};
                zFuncs_mine{i} = zFuncs{2*(i-1)+1}; %OBS could change
        end

        drops.Nps(i,1) = inputstruct.drops{i}.nps;
        drops.lambdas(i,1) = inputstruct.drops{i}.lambda;
        globaldata.lambda(i,1) = inputstruct.drops{i}.lambda;

        surfact.Pe(i,1) = inputstruct.surfact{i}.interf.Pe;
        surfact.E(i,1) = inputstruct.surfact{i}.interf.E;

        globaldata.surfact.rhoinit(i) = inputstruct.surfact{i}.interf.rhoinit;
        globaldata.surfact.E(i) = surfact.E(i,1);
        globaldata.surfact.Pe(i) = surfact.Pe(i);
    end
    zFuncs = zFuncs_mine;


    %The area of the bubbles. Is used as a crude error check.
    [As,origL,~] = get_area(zFuncs,1000);
    origA = sum(As);
    %The original npanel/panel length ratio.
    origq = drops.Nps./origL;
    globaldata.origA = origA;

    globaldata.origq = origq;

    % Discretize the boundary according to the trapezoidal rule with the
    % discretization points being equidistant in arclength.
    drops.idx = zeros(length(zFuncs),2);
    drops.z = zeros(sum(drops.Nps)*drops.pointsPerPanel,1);
    drops.W = zeros(sum(drops.Nps)*drops.pointsPerPanel,1);
    surfact.rho = zeros(sum(drops.Nps)*drops.pointsPerPanel,1);
    surfact.rhok = surfact.rho;
    surfact.rhomass0 = zeros(length(zFuncs),1);
    drops.lambda = zeros(2*sum(drops.Nps)*drops.pointsPerPanel,1);
    drops.bubble = zeros(2*sum(drops.Nps)*drops.pointsPerPanel,1);
    drops.area = origA;
    ptr = 1;

    for j = 1:length(zFuncs)
        [zn,Wn] = Traparcldisc2(@(T) zFuncs_mine{j}(T),drops.pointsPerPanel*drops.Nps(j), ...
            0,drops.pointsPerPanel);
        drops.z(ptr:ptr+drops.pointsPerPanel*drops.Nps(j)-1) = zn;
        surfact.rho(ptr:ptr+drops.pointsPerPanel*drops.Nps(j)-1) = globaldata.surfact.rhoinit(j);
        surfact.rhok(ptr:ptr+drops.pointsPerPanel*drops.Nps(j)-1) =  ...
            1/(drops.pointsPerPanel*drops.Nps(j))*fftshift(...
            fft(surfact.rho(ptr:ptr+drops.pointsPerPanel*drops.Nps(j)-1)));
        equalarc = origL(j)/(2*pi);
        drops.lambda(2*(ptr-1)+(1:2*drops.pointsPerPanel*drops.Nps(j))) = drops.lambdas(j);
        drops.bubble(2*(ptr-1)+(1:2*drops.pointsPerPanel*drops.Nps(j))) = j-1;
        drops.W(ptr:ptr+drops.pointsPerPanel*drops.Nps(j)-1) = Wn(1);
        surfact.rhomass0(j,1) = sum(surfact.rho(ptr:ptr+ ...
            drops.pointsPerPanel*drops.Nps(j)-1).*equalarc...
            .*drops.W(ptr:ptr+drops.pointsPerPanel*drops.Nps(j)-1));
        drops.idx(j,:) = [ptr ptr+drops.pointsPerPanel*drops.Nps(j)-1];
        ptr = ptr + drops.pointsPerPanel*drops.Nps(j);
    end

    surfact.eqstate = globaldata.surfact.eqstate;
    surfact.sigma = compute_surftension(surfact.rho,surfact.E,drops.idx,surfact.eqstate);
    [zp,zpp] = fft_diff(drops.z,drops.idx);
    for j=1:drops.nbrDrops
        surfact.rhocons(j,1) = 0;
        surfact.rhoconsTOT(j,1) = 0;
        surfact.rhomass(j,1) = surfact.rhomass0(j,1);

        I = drops.idx(j,1):drops.idx(j,2);
        surfact.salpha(j,1) = sum(drops.W(I).*abs(zp(I)))/2/pi;
        surfact.kappa(I,1) = imag(zpp(I)./zp(I))./abs(zp(I));

        %         drops.beta(I,1) = (1-drops.lambda(I))./(1+drops.lambda(I));
        %         drops.gamma(I,1) = 1./(1+drops.lambda(I));
    end
    drops.beta = (1-drops.lambda)./(1+drops.lambda);
    drops.gamma = 1./(1+drops.lambda);


end

% Mod for restart files for r1 and r2
if strcmp(ipF, 'indata/restart_r1_Qind50.m')
   drops.beta = (1-drops.lambda)./(1+drops.lambda);
   drops.gamma = 1./(1+drops.lambda);
   surfact.eqstate = 'linear';
end


% Mod to add surfactant concentration to clean drops for swissroll
if strcmp(ipF,'indata/swissroll3_indata.m') || strcmp(ipF,'indata/restart_swissroll3_indata.m')

    initSurf = 0; upSurf = 0;
    switch ipF
        case 'indata/swissroll3_indata.m'
            initSurf = 1;
        case 'indata/restart_swissroll3_indata.m'
            upSurf = 1;
    end

    flds={'z','W','lambda','bubble','beta','gamma'};
    for fld=1:length(flds)
        dropwise.(flds{fld})=cell(drops.nbrDrops,1);
    end

    dropwise.Nps=zeros(drops.nbrDrops,1);
    dropwise.lambdas=zeros(drops.nbrDrops,1);
    dropwise.idx=zeros(drops.nbrDrops,2);

    dropwise_us = dropwise;
    old_drops = drops;

    % Split up
    for drp=1:drops.nbrDrops
        dropwise.z{drp} = drops.z(drops.idx(drp,1):drops.idx(drp,2));
        dropwise.W{drp} = drops.W(drops.idx(drp,1):drops.idx(drp,2));

        dropwise.lambda{drp} = ones(drops.Nps(drp)*drops.pointsPerPanel*2,1)*drops.lambdas(drp);
        dropwise.bubble{drp} = ones(drops.Nps(drp)*drops.pointsPerPanel*2,1);
        dropwise.beta{drp} = (1-dropwise.lambda{drp})./(1+dropwise.lambda{drp});
        dropwise.gamma{drp} = 1./(1+dropwise.lambda{drp});

        dropwise.lambdas(drp) = drops.lambdas(drp);
        dropwise.Nps(drp) = drops.Nps(drp);
        dropwise.idx(drp,:) = drops.idx(drp,:);
    end

    % Upsample and create surfactant concentration
    ptsTot = 0;
    for j=1:drops.nbrDrops
        zj = dropwise.z{j};
        N = length(zj);
        zjk = 1/N*fftshift(fft(zj));
%         zjk(1:floor(N/4)) = 0;
%         zjk(ceil(3*N/4):end) = 0;
        zjk(abs(zjk)<1e-8) = 0;

        upFac = inputstruct.restartupsamplefactor(j);
        newNps = ceil(upFac*N/drops.pointsPerPanel);
        newN = newNps*drops.pointsPerPanel;
        sdiff = newN-length(zjk);

        if sdiff > 0
            new_zjk = [zeros(floor(newN-N)/2,1); zjk; zeros(floor(newN-N)/2,1)];
            zj = newN*ifft(ifftshift(new_zjk));
        end

        dropwise_us.z{j} = zj;

        dropwise_us.Nps(j) = newNps;

        dropwise_us.lambdas(j) = dropwise.lambdas(j);

        if j==1
            dropwise_us.idx(j,1) = 1;
        else
            dropwise_us.idx(j,1) = dropwise_us.idx(j-1,2) + 1;
        end
        dropwise_us.idx(j,2) = dropwise_us.idx(j,1) + newN - 1;

        dropwise_us.lambda{j} = ones(dropwise_us.Nps(j)*drops.pointsPerPanel*2,1)*drops.lambdas(j);
        dropwise_us.bubble{j} = (j-1)*ones(dropwise_us.Nps(j)*drops.pointsPerPanel*2,1);
        dropwise_us.beta{j} = (1-dropwise_us.lambda{j})./(1+dropwise_us.lambda{j});
        dropwise_us.gamma{j} = 1./(1+dropwise_us.lambda{j});

        dropwise_us.W{j} = 2*pi/newN*ones(newN,1);

        ptsTot = ptsTot + newN;
    end

    % Assemble
    for j=1:drops.nbrDrops
        if j==1
            for fld=1:length(flds)
                drops_us.(flds{fld}) = dropwise_us.(flds{fld}){j};
            end
        else
            for fld=1:length(flds)
                drops_us.(flds{fld}) = [drops_us.(flds{fld});dropwise_us.(flds{fld}){j}];
            end
        end
    end
    drops_us.idx=dropwise_us.idx;
    drops_us.Nps=dropwise_us.Nps;
    drops_us.lambdas=dropwise_us.lambdas;
    drops_us.pointsPerPanel = drops.pointsPerPanel;
    drops_us.nbrDrops = drops.nbrDrops;
    drops_us.area = drops.area;
    drops_us.oldarea = drops.oldarea;

    drops = drops_us;


    % Surfactants
    if initSurf || upSurf
        surfact_us=surfact;
        flds={'rho','rhok','sigma','kappa'};
        for fld=1:length(flds)
            surfwise.(flds{fld})=cell(drops.nbrDrops,1);
            surfact_us.(flds{fld})=[];
        end

        for j=1:drops.nbrDrops

            if initSurf
                surfact_us.Pe(j) = inputstruct.surfact{j}.interf.Pe;
                surfact_us.E(j) = inputstruct.surfact{j}.interf.E;
                globaldata.surfact.rhoinit(j) = inputstruct.surfact{j}.interf.rhoinit;

                surfwise.rho{j} = globaldata.surfact.rhoinit(j)*ones(size(dropwise_us.z{j}));
                surfwise.rhok{j} = 1/(drops.Nps(j)*16)*fftshift(fft(surfwise.rho{j}));

            else
                if upSurf
                    surfact_us.Pe(j) = surfact.Pe(j);
                    surfact_us.E(j) = surfact.E(j);

                    rhok = surfact.rhok(old_drops.idx(j,1):old_drops.idx(j,2));
                    upFac = inputstruct.restartupsamplefactor(j);
                    N = old_drops.Nps(j)*drops.pointsPerPanel;
                    newNps = ceil(upFac*N/drops.pointsPerPanel);
                    newN = newNps*drops.pointsPerPanel;
                    sdiff = newN-length(rhok);

                    if sdiff > 0
                        new_rhok = [zeros(floor(newN-N)/2,1); rhok; zeros(floor(newN-N)/2,1)];
                        new_rho = newN*ifft(ifftshift(new_rhok));
                    else
                        new_rhok = rhok;
                        new_rho = surfact.rho(old_drops.idx(j,1):old_drops.idx(j,2));
                    end

                    surfwise.rho{j} = new_rho;
                    surfwise.rhok{j} = new_rhok;

                end
            end

            globaldata.surfact.E(j) = surfact_us.E(j);
            globaldata.surfact.Pe(j) = surfact_us.Pe(j);

            [zp,zpp] = fft_diff(dropwise_us.z{j},[1 drops.Nps(j)*16]);

            surfact_us.salpha(j) = sum(dropwise_us.W{j}.*abs(zp))/2/pi;
            surfact_us.rhomass(j,2) = abs(sum(surfwise.rho{j}.*surfact_us.salpha(j).*dropwise_us.W{j}));
            surfwise.kappa{j} = imag(zpp./zp)./abs(zp);

            surfwise.sigma{j} = compute_surftension(surfwise.rho{j},surfact_us.E(j),[1 drops.Nps(j)*16],surfact_us.eqstate);

        end
        surfact_us.eqstate = globaldata.surfact.eqstate;
        surfact_us.rhomass0 = surfact_us.rhomass(:,2);
        globaldata.surfact.eqstate = inputstruct.surfact{1}.interf.eqstate;

        % Assemble
        for j=1:drops.nbrDrops
            if j==1
                for fld=1:length(flds)
                    surfact_us.(flds{fld}) = surfwise.(flds{fld}){j};
                end
            else
                for fld=1:length(flds)
                    surfact_us.(flds{fld}) = [surfact_us.(flds{fld});surfwise.(flds{fld}){j}];
                end
            end
        end
        surfact = surfact_us;
    end


end



end
