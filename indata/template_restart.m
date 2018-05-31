% Read in and restart. Change final time.

clear all;

inputstruct.restartfile = 'template_outdata_t0p01';

% Time parameters
inputstruct.timedisc.Tint = [0 2]; %time-interval
inputstruct.timedisc.dtmax = 0.01; %largest step size for adaptive time-step
inputstruct.timedisc.Nt = (inputstruct.timedisc.Tint(2)-inputstruct.timedisc.Tint(1))/inputstruct.timedisc.dtmax; %non-adaptive number of adaptive substeps
inputstruct.actionintervals.Ns = 10; %number of steps between save-points
inputstruct.actionintervals.Na = 1; %number of steps between adaptivity checks
inputstruct.actionintervals.Nf = inputstruct.timedisc.Nt; %number of steps between filtering;
inputstruct.actionintervals.Np = 100;

inputstruct.timedisc.hdivisor = 10; %initial estimate of number of steps per substep


% Velocity parameters
load('validsteady_N2400');  
inputstruct.farfield.Q = 0.1; %Far field velocity, irrotational flow =0,0.1
inputstruct.farfield.G = 0; %Far field velocity, shear flow

% Tolerance parameters
inputstruct.steadyTol = 1e-7;
inputstruct.timeTol = 0.5e-6; %In case of adaptive time stepper
inputstruct.circleTol = 0;
inputstruct.doadaptivity = 1;

% Interface parameters
inputstruct.nbrDrops = 1;
inputstruct.drops{1}.shape = 'ellipse';
inputstruct.drops{1}.aspect = 1;
inputstruct.drops{1}.centre = 0;
inputstruct.drops{1}.lambda = 1;
inputstruct.drops{1}.nps = 250;
inputstruct.drops{1}.pointsperpanel = 16;
%... if more drops, add more of these here

% Surfactant parameters
inputstruct.surfact{1}.interf.rhoinit = 1;
inputstruct.surfact{1}.interf.E = 0.5;
inputstruct.surfact{1}.interf.Pe = inf;
inputstruct.surfact{1}.interf.eqstate = 'linear';


try load(['../output/data/' inputstruct.restartfile])
    inputstruct.isrestartfile = 1;
    inputstruct.timedisc.Tint(1) = simdata.t;
    inputstruct.restartfile = ['output/data/' inputstruct.restartfile];
    disp(['Load restartfile from ' inputstruct.restartfile])
catch
    disp(['Restartfile does not exist: ' inputstruct.restartfile])
    inputstruct.isrestartfile = 0;
end


