%Set up standard test bed (1)

% Time parameters
input.timedisc.Tint = [0 1];
input.timedisc.Nt = 10;
input.timedisc.hdivisor = 10;

% Velocity parameters
input.farfield.Q = 0; %Far field velocity, irrotational flow =0,0.1
input.farfield.G = 0; %Far field velocity, shear flow

% Tolerance parameters
input.steadyTol = 1e-8;
input.timeTol = 1e-5; %In case of adaptive time stepper
input.circleTol = 0;
input.doadaptivity = 0;

% Interface parameters
input.nbrDrops = 1;
input.drops{1}.shape = 'ellipse';
input.drops{1}.aspect = 1;
input.drops{1}.centre = 0;
input.drops{1}.lambda = 1;
input.drops{1}.nps = 32;
input.drops{1}.pointsperpanel = 16;
%... if more drops, add more of these here

% Surfactant parameters
input.surfact{1}.interf.rhoinit = 1;
input.surfact{1}.interf.E = 0.5;
input.surfact{1}.interf.Pe = inf;
input.surfact{1}.interf.eqstate = 'linear';



