% InitParamMaster
%
% Parameter inializer for RunHardRod

% Main flagMaster
flagMaster.SaveMe       = 1; % Saving
flagMaster.Verbose      = 0; % Prints run times
flagMaster.AnisoDiff    = 1; % Aniso = 1. Iso = 0
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

particleMaster.lMaj = 1;  % Length along the major axis
particleMaster.lMin = 0;  % Length along the minor axis
particleMaster.vD   = [0]; % Driving velocity (vec)
particleMaster.mob  = 1; % mobility

systemMaster.Nx = [64]; % Gridpoints in x dir (vec)
systemMaster.Ny = [64]; % Gridpoints in y dir (vec) 
systemMaster.Nm = [64]; % Gridpoints in angle (vec)
systemMaster.bc = [1.45]; % Scaled concentration (vec)
systemMaster.Lx = [10];  % Box length (vec)
systemMaster.Ly = [10];  % Box length (vec)
systemMaster.Lphi = 2 * pi;
systemMaster.Tmp = 1;   % Temperature

if flagMaster.AnisoDiff
  particleMaster.mobPar  = 2*particleMaster.mob; 
else
  particleMaster.mobPar = particleMaster.mob; 
end
particleMaster.mobPerp   = particleMaster.mob;
particleMaster.mobRot   = 6 * particleMaster.mob / particleMaster.lMaj^2;

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
%phi
rhoInitMaster.aPhif = 0; % Gauss amp in phi fraction of concentration
rhoInitMaster.varPhi = 0; % Variance of gaussian in phi
%x
rhoInitMaster.aXf = 0; % Gauss amp in x as fraction of concentration
rhoInitMaster.varX  = 0; % Variance of gaussian in x
rhoInitMaster.centerX = 0; % Center of gaussian in x
%y
rhoInitMaster.aYf = 0; % Gauss amp in x as fraction of concentration
rhoInitMaster.varY  = 0; % Variance of gaussian in x
rhoInitMaster.centerY = 0; % Center of gaussian in x
% Calculated stuff- fix times, etc.
% Change odd gridspacings to even unless it's one. 
% L=1 for N=1 is for integrations
if systemMaster.Nx == 1 
  systemMaster.Lx = 1;
  rhoInitMaster.NumModesX = 0;
else
  systemMaster.Nx = systemMaster.Nx + mod( systemMaster.Nx, 2 );
end
if systemMaster.Ny == 1 
  systemMaster.Ly = 1;
  rhoInitMaster.NumModesY = 0;
else
  systemMaster.Ny = systemMaster.Ny + mod( systemMaster.Ny, 2 );
end
if systemMaster.Nm == 1 
  systemMaster.Lphi = 1;
  rhoInitMaster.NumModesM = 0;
else
  systemMaster.Nm = systemMaster.Nm + mod( systemMaster.Nm, 2 );
end

% For now, fix flags and turn off traditional analysis if Nm = 1
if systemMaster.Nm == 1
  flagMaster.MakeOP = 0;
  flagMaster.MakeMovies = 0;
  particleMaster.mobRot = 0;
  particleMaster.mobPar = particleMaster.mob;
  particleMaster.mobPerp = particleMaster.mob;
  rhoInitMaster.IntCond = 0;
end

% Don't perturb more more than you are allowed to
if( rhoInitMaster.NumModesX >= systemMaster.Nx / 2 )
  rhoInitMaster.NumModesX = floor(systemMaster.Nx / 2) - 1; 
end
if( rhoInitMaster.NumModesY >= systemMaster.Ny / 2 )
  rhoInitMaster.NumModesY = floor(systemMaster.Ny / 2) - 1;
end
if( rhoInitMaster.NumModesM >= systemMaster.Nm / 2 )
  rhoInitMaster.NumModesM = floor(systemMaster.Nm / 2) - 1; 
end

%  Make sure variance isn't zero if doing polar
if rhoInitMaster.varX ~= 0
  if varX == 0
    rhoInitMaster.varX = systemMaster.Lx/2; 
  end
end
if rhoInitMaster.varY ~= 0
  if varY == 0
    rhoInitMaster.varY = systemMaster.Ly/2; 
  end
end
if rhoInitMaster.varPhi ~= 0
  if varPhi == 0
    rhoInitMaster.varPhi = systemMaster.Lphi/2; 
  end
end
% Scale ss_epsilon by delta_t. Equivalent to checking d rho /dt has reached
% steady state instead of d rho
timeMaster.ss_epsilon_dt = timeMaster.ss_epsilon .* timeMaster.dt;
% Make sure step method 6 isn't on ani = 0
if flagMaster.AnisoDiff == 0 && flagMaster.StepMeth == 6
  flagMaster.StepMeth = 0;
end
% Fix Ls if we want the box to be square
if flagMaster.SquareBox == 1
  systemMaster.L_box = unique( [systemMaster.Lx systemMaster.Ly] );
  systemMaster.Lx = systemMaster.L_box; 
  systemMaster.Ly = systemMaster.L_box; 
end
% Fix Lx is we want all Ns to be the same
if flagMaster.AllNsSame == 1
  Nvec = unique( [systemMaster.Nx systemMaster.Ny systemMaster.Nm] );
  systemMaster.Nx = Nvec;
  systemMaster.Ny = Nvec; 
  systemMaster.Nm = Nvec;
end
% Concentration and rod stuff
particleMaster.b  = particleMaster.lMaj^2/pi;               % Average excluded volume per particle
particleMaster.mobPos = particleMaster.mob;
% Make OP if making movies 
if flagMaster.MakeMovies == 1; Flag.MakeOP = 1; end % if make movie, make OP first
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
