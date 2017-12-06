% InitParamMaster
% Parameter inializer for runHardRod
% Main flagMaster
flagMaster.SaveMe       = 1; % Saving
flagMaster.parforFlag   = 0; % Saving
flagMaster.Verbose      = 0; % Prints run times
flagMaster.DiagLop      = 1; % Diag operator = 1. off diag = 0
flagMaster.MakeOP       = 1; % No OPs if save is zero
flagMaster.AllNsSame    = 0; % Sets all gridptns to be the same
flagMaster.SquareBox    = 0; % Forces box to be square
flagMaster.StepMeth     = [6]; % Stepping (integrating) method (vec)
flagMaster.rndStrtUpSeed = 1; % start with a random seed (1) or startup seed (0)
flagMaster.scaleDt = 1; % scale dt by velocity
% 0: AB1 1: AB2 2: HAB1 3: HAB2 4: BHAB1 5: BHAB2 6: phiV- Aniso EE-Iso
% movie flags
movieFlagMaster.analysis = 0; % Turn off all plots
movieFlagMaster.makeMovie = 0; % movies
movieFlagMaster.plotInset = 0; % plot inset with movies
movieFlagMaster.plotMax = 0; % max vs time
movieFlagMaster.plotAmp = 0; % plot Ft amplitudes
movieFlagMaster.plotSlice = 0; % plot crysal structure
movieFlagMaster.plotCrystal = 0; % plot crysal structure
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
particleMaster.nlDr = [0]; % off if 0, if ~0, rhoMax value
particleMaster.nlDiff = {0, []}; % off if 0, if ~0, rhoMax value
%  Long range interaction type in MF: 
% {'ss', 'mf/vir', es1, es2, lr1, lr2}  (softshoulder2d), 
% {'pa', 'mf/vir', es1} (polaralign2d), 
% {'pag', 'mf/vir', es1, ls1} (polaraligngauss2d), 
% {'de', 'mf/vir', es1, ls1} (decay exponential)
% {'gauss', 'mf/vir', es1, ls1} (gaussian)
% different cells are different runs
particleMaster.interactLrV = { }; 
% External potential (cell of cells): 
%  { {'linV', dim, a}, {'quadV', dim, k}, {'nemV', es1, phase} }
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
systemMaster.noise = {0, 0}; % Noise strength {pos, rot}, scaled by gridspacing
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
% if not specified, input should be a double
% {'iso'}: isotropic
% {'nem', shiftAngle}: nematic, nematis
% {'eq'}: equilibrium (hardrod, softshoulder)
% {'load', loadStr, loadpath}: load density default path: ./src/InitialDensities/SavedRhos/
% {'delP', dirAngle}: delta function in polar order
% {'crys', latticeSpacing}: crystal
% {'gauss', amp, var1, center1,  var2, center2, var3, center3} 
% {'lorenz', amp, width1, center1,  width2, center2, width3, center3}
rhoInitMaster.intCond = {'iso'};
% Weight of the spatial sinusoidal perturbation. %
% Perturbation weight is a fraction of the isotropic density
% If about a nematic, code will correct for negative densities.
% perturb = 
% {'none'}: no perturbations
% {'pw', int randomAmp, weight, int numModes1, int numModes2, int numModes3 }: planewave perturbations.
% {'gauss', amp1, var1, center1,  var2, center2, var3, center3}: gaussian perturbation. homoFlag for concentration
% {'lorenz', amp, width1, center1,  width2, center2, width3, center3}
rhoInitMaster.perturb = { {'pw', 1, 1e-3, 8, 8, 8} };
% Save the Params

save('Params.mat','particleMaster','systemMaster',...
  'runMaster','timeMaster','flagMaster','rhoInitMaster')
