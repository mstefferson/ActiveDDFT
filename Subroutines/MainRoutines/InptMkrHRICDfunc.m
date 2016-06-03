% Input maker. Takes concentration and velocity as inputs

function [FileInpt] = InptMkrHRICDfunc(bc,v0,L_box,N,t_tot,trial)
% Input creater for HR2DrotMainDrIDCube

%%%%%% Turn on/off interactions%%%%%%%%%
Interactions = 1; 
SaveMe       = 0;   
MakeMovies   = 0; % Movies won't run if save is zero
MakeOP       = 0;

%%%%%%%%%%%%% Box and Rod Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx      = N(1);Ny = N(2); Nm = N(3);
L_rod   = 1;                  % Length of the rods
Lx      = L_box;               % Box length
Ly      = L_box;               % Box length
AspctRt = 8;                  % L / W
% v0      = 0;                  %Driving velocity

%%%%%%%%%%%%%%%Time recording %%%%%%%%%%%%%%%%%%%%%%%%%%
dt     = 1e-3; %time step
% t_tot       = 10;   %total time
t_record    = t_tot/100; %time interval for recording dynamics
ss_epsilon  = 1e-10;                          %steady state condition

% The number of k-modes above and below k = 0 added as a perturbation 
NumModesX   = 6;
NumModesY   = 6;
NumModesM   = 6;


%%%%%%%%%%%%%%%%%%%%% Physical Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmp      = 1;            % Temperature
% mobility
Mob_same = 1;
Mob_pos  = Mob_same; 
Mob_rot  = Mob_same;

%%%%%%%%% Initial density parameters%%%%%%%%%%%%%%%%%%
 % Weight of the spatial sinusoidal perturbation. %
 % Perturbations added to rho(i,j,k) = 1. Must be small
WeightPos      = 1e-2;   
WeightAng      = 1e-2;    
Random         = 1;       % Random perturbation coeffs
% bc             = 1.6;    % Dimensionless  scaled concentration

% Type of initial Condition
IntGauss   = 0;
IntPw      = 0;
IntSepPw   = 0;
IntEqPw    = 0;    % Distribution from perturbing equil. dist.
IntEqSepPw = 0;
IntLoad    = 0;
IntNemPw   = 1;    % Distribution from perturbing equil. dist.

% Save a string saying what you want
[IntDenType, IntDenIndicator] = ...
         IntDenIndicatorMaker(IntGauss,IntPw,IntSepPw,IntEqPw,...
         IntEqSepPw,IntLoad,IntNemPw);
     
if IntLoad
    IntDenName = sprintf('DenBlow2');
else
    IntDenName = sprintf('Irrelevant');
end

% Concentration and rod stuff
b       = L_rod^2/pi;               % Average excluded volume per particle
c       = bc / b;                   % Concentration
Norm    = c * Lx * Ly;              % number of particles
Diam    = L_rod / AspctRt;              % "Diameter" aka width of the rods

% Turn movies off is Save is off
if SaveMe == 0; MakeMovies = 0;end
if MakeMovies == 1; MakeOP = 1; end % if make movie, make OP first

% Turn off Drive parameter  is v0 = 0;
if v0  == 0; Drive = 0; else Drive = 1;end

% Parameter vectors
Paramtmp = [trial Interactions Drive MakeOP MakeMovies SaveMe Nx Ny Nm...
    Lx Ly L_rod Tmp Norm WeightPos WeightAng Random ...
    NumModesX NumModesY NumModesM bc Mob_pos Mob_rot v0];
Timetmp  = [dt t_record t_tot ss_epsilon];

% Make the output directory string and input file
FileDir = ...
    sprintf('HRdiffIDC_N%i%i%i_bc%.1f_Int%i_v%.1f_%st%i',...
    Nx,Ny,Nm,bc,Interactions,v0,IntDenIndicator,trial);
FileInpt = ...
    sprintf('Inpt_N%i%i%i_bc%.1f_Int%i_v%.1f_%st%i.txt', ...
    Nx,Ny,Nm,bc,Interactions,v0,IntDenIndicator,...
    trial);

Where2SavePath    = sprintf('%s/%s/%s',pwd,'Outputs',FileDir);

fid = fopen(FileInpt,'wt');  % Note the 'wt' for writing in text mode
fprintf(fid,'%s\n%s\n%s\n',FileDir,Where2SavePath,IntDenType);
fprintf(fid,'Param\t');
fprintf(fid,'%e\t',Paramtmp);
fprintf(fid,'\n');
fprintf(fid,'Time\t');
fprintf(fid,'%e\t',Timetmp);

fclose(fid);






