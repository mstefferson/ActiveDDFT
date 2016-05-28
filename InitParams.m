% Anisotropic diffusion
% Input creater for HR2DrotMainDr

cd ~/DDFT/HardRodML
CurrentDir = pwd;
addpath( genpath( CurrentDir ) );

% Main Flags
Flags.SaveMe       = 1; % Saving
Flags.AnisoDiff    = 1; % Aniso = 1. Iso = 0
Flags.Interactions = 1; % Int = 1. No Int = 0
Flags.MakeMovies   = 1; % No movies if save is zero
Flags.MakeOP       = 1; % No OPs if save is zero
% Stepping---0: AB1 1: AB2 2: HAB1 3: HAB2 4: BHAB1 5: BHAB2 6: phiV- Aniso EE-Iso
Flags.StepMeth     = 1;

%%%%%% Parameters %%%%%%%%%%%%%%%%%
% Can be vector if (vec) is in comment
ParamMaster.trial   = 200;  % Trial Indicator
ParamMaster.Nx      = 64; % Gridpoints in x dir (vec)
ParamMaster.Ny      = 64; % Gridpoints in y dir (vec) 
ParamMaster.Nm      = 64; % Gridpoints in angle (vec)
ParamMaster.bc      = 1.45; % Scaled concentration (vec)
ParamMaster.L_rod   = 1;  % Length of the rods
ParamMaster.Lx      = 10*L_rod;  % Box length (vec)
ParamMaster.Ly      = 10*L_rod;  % Box length (vec)
ParamMaster.vD      = 30.0; % Driving velocity (vec)
ParamMaster.Tmp      = 1;            % Temperature
% mobility
Mob = 1;
if Flags.Aniso; ParamMaster.Mob_par  = 2*Mob; else; ParamMaster.Mob_par = Mob; end;
ParamMaster.Mob_perp   = Mob;
ParamMaster.Mob_rot   = 6 * Mob / L_rod^2;

%%%%%%%%%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%%
TimeObj.delta_t     = 1e-3; % time step
TimeObj.t_record    = 0.1; % time interval for recording dynamics
TimeObj.t_tot       = 10.0; % total time
TimeObj.ss_epsilon  = 1e-8;  % steady state condition

%%%%%%%%% Initial Condition %%%%%%%%%%%%%%%%%%%%%
RhoInit.IntCond     = 1;
RhoInit.NumModesX   = 8;
RhoInit.NumModesY   = 8;
RhoInit.NumModesM   = 8;
% Don't perturb more more than you are allowed to
if( NumModesX >= Nx / 2 ); NumModesX = floor(Nx / 2) - 2; end;
if( NumModesY >= Ny / 2 ); NumModesY = floor(Ny / 2) - 2; end;
if( NumModesM >= Nm / 2 ); NumModesM = floor(Nm / 2) - 2; end;

% Weight of the spatial sinusoidal perturbation. %
% Perturbations added to rho(i,j,k) = 1. Must be small
RhoInit.WeightPos   = 1e-3;
RhoInit.WeightAng   = 1e-3;
RhoInit.RandomAmp   = 1;       % Random perturbation coeffs

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

% Concentration and rod stuff
b       = ParamMaster.L_rod^2/pi;               % Average excluded volume per particle
ParamMaster.c       = bc / b;                   % Concentration
%%%%%%%%%%%%%% DO THIS IN RUN %%%%%%%%%%%%%
ParamMaster.Norm    = c * Lx * Ly;              % number of particles

% Fix the time
[TimeObj.t_tot,TimeObj.N_time,TimeObj.t_rec,TimeObj.N_rec,TimeObj.N_count]= ...
        TimeStepRecMaker(delta_t,t_tot,t_record);

% Turn movies off is Save is off
if Flags.SaveMe == 0; Flags.MakeOP = 0; Flags.MakeMovies = 0;;end
if Flags.MakeMovies == 1; Flag.MakeOP = 1; end % if make movie, make OP first
if ParamMaster.vD  == 0; Flags.Drive = 0; else Flags.Drive = 1;end

% Save the Params
save('Params.mat','ParamMaster','TimeObj','Flags','RhoInit')
