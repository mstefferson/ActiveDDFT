% InitParamMaster
% Parameter inializer for runHardRod
% Main flagMaster
flagMaster.SaveMe       = 1; % Saving
flagMaster.Verbose      = 0; % Prints run times
flagMaster.DiagLop      = 1; % Diag operator = 1. off diag = 0
flagMaster.MakeMovies   = 0; % No movies if save is zero
flagMaster.MakeOP       = 1; % No OPs if save is zero
flagMaster.AllNsSame    = 0; % Sets all gridptns to be the same
flagMaster.SquareBox    = 0; % Forces box to be square
flagMaster.StepMeth     = [6]; % Stepping (integrating) method (vec)
flagMaster.rndStrtUpSeed = 1; % start with a random seed (1) or startup seed (0)
flagMaster.scaleDt = 1; % scale dt by velocity
% 0: AB1 1: AB2 2: HAB1 3: HAB2 4: BHAB1 5: BHAB2 6: phiV- Aniso EE-Iso
%%%%%% Parameters %%%%%%%%%%%%%%%%%
% Can be vector if (vec) is in comment
runMaster.numTrial =  1; % number of trials
runMaster.trialID = 1; % Trial Indicator
runMaster.runID   = 1; % Starting Run indicator *_trial.runID
%%%%%%%%%% Particles %%%%%%%%%%
particleMaster.numTypes = 1; % number of types
particleMaster.lMaj = 1;  % Length along the major axis
particleMaster.lMin = 0;  % Length along the minor axis
particleMaster.vD   = [0]; % Driving velocity (vec)
particleMaster.mob  = 1; % mobility
particleMaster.type = 'rods'; % rods, disks, spheres 
particleMaster.interHb = 'mayer'; % Hard body interactions type mayer, spt
%  Long range interaction type in MF: 
% {'ss', 'mf/vir', es1, es2, lr1, lr2}  (softshoulder2d), 
% {'pa', 'mf/vir', es1} (polaralign2d), 
% {'pag', 'mf/vir', es1, ls1} (polaraligngauss2d), 
% {'de', 'mf/vir', es1, ls1} (decay exponential)
% {'gauss', 'mf/vir', es1, ls1} (gaussian)
% different cells are different runs
particleMaster.interactLrV = { }; 
% External potential (cell of cells): 
%  { {'linV', dim, a}, {'quadV', dim, a}, {'nemV', es1, phase} }
particleMaster.externalV = { }; % External potential 
%%%%%%%%% System %%%%%%%%%%%%%%%%
systemMaster.n1 = [64]; % Gridpoints in x dir (vec)
systemMaster.n2 = [64]; % Gridpoints in y dir (vec) 
systemMaster.n3 = [64]; % Gridpoints in angle (vec)
systemMaster.bc = [1.45]; % Scaled concentration (vec)
systemMaster.l1 = [10];  % Box length (vec)
systemMaster.l2 = [10];  % Box length (vec)
systemMaster.l3 = 2 * pi;
systemMaster.tmp = 1;   % Temperature
%%%%%%%%%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%%
timeMaster.dt         = 1e-3; % time step
timeMaster.t_rec      = 0.1;  % time elapsed before recording
timeMaster.t_write    = 0.2;  % time elapsed before writing to file
timeMaster.t_tot      = 1.0;  % total run time
timeMaster.ss_epsilon = 1 * 10^(-5); % steady state condition dRho
timeMaster.scaleDt = flagMaster.scaleDt;
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
% 7: Crystal
% 8: Lorenzian
rhoInitMaster.RandomAmp = 1; % Random perturbation coeffs
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
% Crystal crystalLattice = [a sigma]. A: 2.26, B:1.21
rhoInitMaster.crystalLattice = [1.21 systemMaster.l1 / 100]; 
% Save the Params
save('Params.mat','particleMaster','systemMaster',...
  'runMaster','timeMaster','flagMaster','rhoInitMaster')
