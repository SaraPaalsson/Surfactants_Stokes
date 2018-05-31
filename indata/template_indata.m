
% List of parameters, defined for stokes_surf
% Template updated: 2016-11-28

clearvars -except ipF

% Time parameters
inputstruct.timedisc.Tint = [0 0.02]; %time-interval
inputstruct.timedisc.dtmax = 0.01; %largest step size for adaptive time-step
inputstruct.timedisc.dt = 0.01;
inputstruct.timedisc.Nt = (inputstruct.timedisc.Tint(2)-inputstruct.timedisc.Tint(1))/inputstruct.timedisc.dt; %non-adaptive number of adaptive substeps
inputstruct.actionintervals.Ns = 10; %number of steps between save-points
inputstruct.actionintervals.Na = 1; %number of steps between adaptivity checks
inputstruct.actionintervals.Nf = inputstruct.timedisc.Nt; %number of steps between filtering;
inputstruct.actionintervals.Np = 1;

inputstruct.timedisc.hdivisor = 1; %initial estimate of number of steps per substep


% Velocity parameters
inputstruct.farfield.Q = 0; %Far field velocity, irrotational flow =0,0.1
inputstruct.farfield.G = 0; %Far field velocity, shear flow

% Tolerance parameters
inputstruct.steadyTol = 0;
inputstruct.timeTol = 0.5e-6; %In case of adaptive time stepper
inputstruct.circleTol = 0;
inputstruct.doadaptivity = 1;

% Interface parameters
inputstruct.nbrDrops = 2;
inputstruct.drops{1}.shape = 'ellipse';
inputstruct.drops{1}.aspect = 2;
inputstruct.drops{1}.centre = 0.72i;
inputstruct.drops{1}.lambda = 0;
inputstruct.drops{1}.nps = 25;
inputstruct.drops{1}.pointsperpanel = 16;
%... if more drops, add more of these here

% Surfactant parameters
inputstruct.surfact{1}.interf.rhoinit = 0;
inputstruct.surfact{1}.interf.E = 0;
inputstruct.surfact{1}.interf.Pe = inf;
inputstruct.surfact{1}.interf.eqstate = 'linear';

% Interface parameters
inputstruct.drops{2}.shape = 'ellipse';
inputstruct.drops{2}.aspect = 2;
inputstruct.drops{2}.centre = -inputstruct.drops{1}.centre;
inputstruct.drops{2}.lambda = 0;
inputstruct.drops{2}.nps = inputstruct.drops{1}.nps;
inputstruct.drops{2}.pointsperpanel = 16;
%... if more drops, add more of these here

% Surfactant parameters
inputstruct.surfact{2}.interf.rhoinit = 0;
inputstruct.surfact{2}.interf.E = 0;
inputstruct.surfact{2}.interf.Pe = inf;
inputstruct.surfact{2}.interf.eqstate = 'linear';

