% InitParamMaster
% Parameter inializer for runHardRod
% Main flagMaster
flagMaster.SaveMe       = 1; % Saving
flagMaster.Verbose      = 0; % Prints run times
flagMaster.DiagLop      = 0; % Diag operator = 1. off diag = 0
flagMaster.Interactions = 1; % Int = 1. No Int = 0
flagMaster.MakeMovies   = 0; % No movies if save is zero
flagMaster.MakeOP       = 1; % No OPs if save is zero
flagMaster.AllNsSame    = 0; % Sets all gridptns to be the same
flagMaster.SquareBox    = 0; % Forces box to be square
flagMaster.StepMeth     = [6]; % Stepping (integrating) method (vec)
flagMaster.rndStrtUpSeed = 1; % start with a random seed (1) or startup seed (0)
% 0: AB1 1: AB2 2: HAB1 3: HAB2 4: BHAB1 5: BHAB2 6: phiV- Aniso EE-Iso
%%%%%% Parameters %%%%%%%%%%%%%%%%%
% Can be vector if (vec) is in comment
runMaster.num_trial =  1; % number of trials
runMaster.trialID = 1; % Trial Indicator
runMaster.runID   = 1; % Starting Run indicator *_trial.runID
%%%%%%%%%% Particles %%%%%%%%%%
particleMaster.numTypes = 1; % number of types
particleMaster.type = 'rods'; % rods, disks, spheres 
particleMaster.interHb = 'mayer'; % Hard body interactions type mayer, spt, fmt
particleMaster.interLr = []; %  Long range interaction type in MF.
particleMaster.interLrEnSc = []; %  Long range interaction energy scale in MF.
particleMaster.externalPot = []; % External potential 
particleMaster.lMaj = 1;  % Length along the major axis
particleMaster.lMin = 0;  % Length along the minor axis
particleMaster.vD   = [0]; % Driving velocity (vec)
particleMaster.mob  = 1; % mobility
%%%%%%%%% System %%%%%%%%%%%%%%%%
systemMaster.n1 = [64]; % Gridpoints in x dir (vec)
systemMaster.n2 = [64]; % Gridpoints in y dir (vec) 
systemMaster.n3 = [64]; % Gridpoints in angle (vec)
systemMaster.bc = [1.45]; % Scaled concentration (vec)
systemMaster.l1 = [10];  % Box length (vec)
systemMaster.l2 = [10];  % Box length (vec)
systemMaster.l3 = 2 * pi;
systemMaster.Tmp = 1;   % Temperature
%%%%%%%%%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%%
timeMaster.dt         = 1e-3; % time step
timeMaster.t_rec      = 0.1;  % time elapsed before recording
timeMaster.t_write    = 0.2;  % time elapsed before writing to file
timeMaster.t_tot      = 1.0;  % total run time
timeMaster.ss_epsilon = 1 * 10^(-5); % steady state condition dRho
%%%%%%%%% Initial Condition %%%%%%%%%%%%%%%%%%%%%
% Key
% The number of k-modes above and below k = 0 added as a perturbation
% Type of Inital Condition
% 0: Plane wave perturbation over an isotropic distribution
% 1: Plane wave perturbation over the equilibrium distribution
% 2: Plane wave perturbation over a nematic distribution
% 3: Load an equilbrium distribution and fit it to your box
% 4: Gaussian perturbations with homogenous concentration 
% 5: Gaussian perturbations with inhomogenous concentration 
% 6: Delta function in polar order 
rhoInitMaster.IntCond   = [1]; % IC indicator (vec)
rhoInitMaster.LoadName  = ''; % File name to load if IC 3
% Weight of the spatial sinusoidal perturbation. %
% Perturbation weight is a fraction of the isotropic density
% If about a nematic, code will correct for negative densities.
rhoInitMaster.RandomAmp = 1; % Random perturbation coeffs
rhoInitMaster.NumModesX = 8; % Perturb # modes x
rhoInitMaster.NumModesY = 8; % Perturb # modes y
rhoInitMaster.NumModesM = 8; % Perturb # modes m
rhoInitMaster.WeightPert = 1e-3; % Weight of perturbations
rhoInitMaster.shiftAngle = 0; % If perturbing about nematic, shift distribution by this much
% Gaussian perturbation 
% rhoInitMaster.gP = [ aYf, varY, centerY, aYf, varY, centerY, aPhif, varPhi, centerPhi ]
rhoInitMaster.gP = [0, systemMaster.l1/2, 0, 0, systemMaster.l1/2, 0, 0, systemMaster.l3/2, 0]; 
% Calculated stuff- fix times, etc.
% Change odd gridspacings to even unless it's one. 
% L=1 for N=1 is for integrations
if systemMaster.n1 == 1 
  systemMaster.l1 = 1;
  rhoInitMaster.NumModesX = 0;
else
  systemMaster.n1 = systemMaster.n1 + mod( systemMaster.n1, 2 );
end
if systemMaster.n2 == 1 
  systemMaster.l2 = 1;
  rhoInitMaster.NumModesY = 0;
else
  systemMaster.n2 = systemMaster.n2 + mod( systemMaster.n2, 2 );
end
if systemMaster.n3 == 1 
  systemMaster.l3 = 1;
  rhoInitMaster.NumModesM = 0;
else
  systemMaster.n3 = systemMaster.n3 + mod( systemMaster.n3, 2 );
end
% Get particle mobility
particleMaster = particleMobInit( particleMaster, flagMaster.DiagLop);
% For now, fix flags and turn off traditional analysis if n3 = 1
if systemMaster.n3 == 1
  flagMaster.MakeOP = 0;
  flagMaster.MakeMovies = 0;
  particleMaster.mobRot = 0;
  particleMaster.mobPar = particleMaster.mob;
  particleMaster.mobPerp = particleMaster.mob;
  rhoInitMaster.IntCond = 0;
end
% Don't perturb more more than you are allowed to
if( rhoInitMaster.NumModesX >= systemMaster.n1 / 2 )
  rhoInitMaster.NumModesX = floor(systemMaster.n1 / 2) - 1; 
end
if( rhoInitMaster.NumModesY >= systemMaster.n2 / 2 )
  rhoInitMaster.NumModesY = floor(systemMaster.n2 / 2) - 1;
end
if( rhoInitMaster.NumModesM >= systemMaster.n3 / 2 )
  rhoInitMaster.NumModesM = floor(systemMaster.n3 / 2) - 1; 
end
%  Make sure variance isn't zero if doing polar
if rhoInitMaster.gP(2) == 0
  rhoInitMaster.gP(2) = systemMaster.l1/2; 
end
if rhoInitMaster.gP(5) == 0
  rhoInitMaster.gP(2) = systemMaster.l1/2; 
end
if rhoInitMaster.gP(8) == 0
  rhoInitMaster.gP(8) = systemMaster.l1/2; 
end
% Scale ss_epsilon by delta_t. Equivalent to checking d rho /dt has reached
% steady state instead of d rho
timeMaster.ss_epsilon_dt = timeMaster.ss_epsilon .* timeMaster.dt;
% Fix Ls if we want the box to be square
if flagMaster.SquareBox == 1
  systemMaster.L_box = unique( [systemMaster.l1 systemMaster.l2] );
  systemMaster.l1 = systemMaster.L_box; 
  systemMaster.l2 = systemMaster.L_box; 
end
% Fix l1 is we want all Ns to be the same
if flagMaster.AllNsSame == 1
  Nvec = unique( [systemMaster.n1 systemMaster.n2 systemMaster.n3] );
  systemMaster.n1 = Nvec;
  systemMaster.n2 = Nvec; 
  systemMaster.n3 = Nvec;
end
% Make OP if making movies 
if flagMaster.MakeMovies == 1; Flag.MakeOP = 1; end % if make movie, make OP first
%Currently, you must save
if flagMaster.MakeOP && flagMaster.SaveMe == 0
  fprintf('Turning on saving, you must be saving to make OPs (due to matfile)\n');
  flagMaster.SaveMe = 1;
end
if particleMaster.vD  == 0; flagMaster.Drive = 0; else flagMaster.Drive = 1;end
% Give time warnings
if timeMaster.dt > timeMaster.t_rec
  fprintf('Recorded interval is shorter then timestep fix before proceeding\n');
end
if timeMaster.t_rec > timeMaster.t_write
  fprintf('File write interval is shorter than recordord interval. Fix it \n');
end
if timeMaster.t_write > timeMaster.t_tot
  fprintf('Total time interval is shorter then record interval\n');
end
% Save the Params
save('Params.mat','particleMaster','systemMaster',...
  'runMaster','timeMaster','flagMaster','rhoInitMaster')
