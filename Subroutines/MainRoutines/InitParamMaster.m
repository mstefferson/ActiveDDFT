% InitParamMaster
%
% Parameter inializer for RunHardRod

%cd ~/DDFT/HardRodML
CurrentDir = pwd;
addpath( genpath( CurrentDir ) );

% Main FlagMaster
FlagMaster.SaveMe       = 1; % Saving
FlagMaster.AnisoDiff    = 1; % Aniso = 1. Iso = 0
FlagMaster.Interactions = 1; % Int = 1. No Int = 0
FlagMaster.MakeMovies   = 1; % No movies if save is zero
FlagMaster.MakeOP       = 1; % No OPs if save is zero
FlagMaster.AllNsSame    = 0; % Sets all gridptns to be the same
FlagMaster.SquareBox    = 0; % Forces box to be square
FlagMaster.StepMeth     = [0]; % Stepping (integrating) method (vec)
% 0: AB1 1: AB2 2: HAB1 3: HAB2 4: BHAB1 5: BHAB2 6: phiV- Aniso EE-Iso

%%%%%% Parameters %%%%%%%%%%%%%%%%%
% Can be vector if (vec) is in comment
ParamMaster.trial   = 1;  % Trial Indicator
ParamMaster.runID   = 1;    % Run indicator *_trial.runID
ParamMaster.Nx      = [64]; % Gridpoints in x dir (vec)
ParamMaster.Ny      = [64]; % Gridpoints in y dir (vec) 
ParamMaster.Nm      = [64]; % Gridpoints in angle (vec)
ParamMaster.bc      = 1.45; % Scaled concentration (vec)
ParamMaster.L_rod   = 1;  % Length of the rods
ParamMaster.Lx      = [10]*ParamMaster.L_rod;  % Box length (vec)
ParamMaster.Ly      = [10]*ParamMaster.L_rod;  % Box length (vec)
ParamMaster.vD      = [0]; % Driving velocity (vec)
ParamMaster.Tmp      = 1;            % Temperature
% mobility
Mob = 1;
if FlagMaster.AnisoDiff; 
  ParamMaster.Mob_par  = 2*Mob; 
else
  ParamMaster.Mob_par = Mob; 
end
ParamMaster.Mob_perp   = Mob;
ParamMaster.Mob_rot   = 6 * Mob / ParamMaster.L_rod^2;

%%%%%%%%%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%%
TimeObj.delta_t     = 1e-3; % time step
TimeObj.t_record    = 0.1; % time interval for recording dynamics
TimeObj.t_tot       = 1.0; % total time
TimeObj.ss_epsilon  = 1e-8;  % steady state condition

%%%%%%%%% Initial Condition %%%%%%%%%%%%%%%%%%%%%
RhoInitMaster.IntCond     = [1]; % IC indicator (vec)
RhoInitMaster.NumModesX   = 8; % Perturb # modes x
RhoInitMaster.NumModesY   = 8; % Perturb # modes y
RhoInitMaster.NumModesM   = 8; % Perturb # modes m
RhoInitMaster.LoadName    = ''; % File name to load if IC 3
% Don't perturb more more than you are allowed to
if( RhoInitMaster.NumModesX >= ParamMaster.Nx / 2 ); 
  RhoInitMaster.NumModesX = floor(ParamMaster.Nx / 2) - 2; 
end;
if( RhoInitMaster.NumModesY >= ParamMaster.Ny / 2 ); 
  RhoInitMaster.NumModesY = floor(ParamMaster.Ny / 2) - 2;
end;
if( RhoInitMaster.NumModesM >= ParamMaster.Nm / 2 );
  RhoInitMaster.NumModesM = floor(ParamMaster.Nm / 2) - 2; 
end;

% Weight of the spatial sinusoidal perturbation. %
% Perturbations added to rho(i,j,k) = 1. Must be small
RhoInitMaster.WeightPos   = 1e-3;
RhoInitMaster.WeightAng   = 1e-3;
RhoInitMaster.RandomAmp   = 1;       % Random perturbation coeffs

% Key
% The number of k-modes above and below k = 0 added as a perturbation
% Type of Inital Condition
% 0: Plane wave perturbation over an isotropic distribution
% 1: Plane wave perturbation over the equilibrium distribution
% 2: Plane wave perturbation over a nematic distribution
% 3: Load an equilbrium distribution
% 4: Seperate plane waves over an isotropic distribution (non-sensical)
% 5: Seperate plane waves over an nematic distribution (non-sensical)
% 6: A gaussian initial condition

% Calculated stuff- fix times, etc.

if FlagMaster.SquareBox == 1; 
  ParamMaster.L_box = unique( [ParamMaster.Lx ParamMaster.Ly] );
  ParamMaster.Lx = ParamMaster.L_box; 
  ParamMaster.Ly = ParamMaster.L_box; 
end

if FlagMaster.AllNsSame == 1; 
  Nvec = unique( [ParamMaster.Nx ParamMaster.Ny ParamMaster.Nm] );
  ParamMaster.Nx = Nvec;
  ParamMaster.Ny = Nvec; 
  ParamMaster.Nm = Nvec;
end


% Concentration and rod stuff
ParamMaster.b  = ParamMaster.L_rod^2/pi;               % Average excluded volume per particle
ParamMaster.Mob_pos = Mob;

% Fix the time
[TimeObj.t_tot,TimeObj.N_time,TimeObj.t_rec,TimeObj.N_rec,TimeObj.N_count]= ...
        TimeStepRecMaker(TimeObj.delta_t,TimeObj.t_tot,TimeObj.t_record);

% Turn movies off is Save is off
if FlagMaster.SaveMe == 0; FlagMaster.MakeOP = 0; FlagMaster.MakeMovies = 0;end
if FlagMaster.MakeMovies == 1; Flag.MakeOP = 1; end % if make movie, make OP first
if ParamMaster.vD  == 0; FlagMaster.Drive = 0; else FlagMaster.Drive = 1;end

% Save the Params
save('Params.mat','ParamMaster','TimeObj','FlagMaster','RhoInitMaster')
