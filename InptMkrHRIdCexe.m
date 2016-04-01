% Input creater for HR2DrotMainIdC
% Isotropic diffusion
cd ~/DDFT/HardRodML 
CurrentDir = pwd;
addpath( genpath( CurrentDir) );

% Now can change number of grid points in the x, y, phi direction
Run  = 1; % Run main from here
Move = 0; % Move files to a nice location

%%%%%%%% Trial %%%%%%%%%%%%
trial    = 10;

% Date
DateTimeStart =  datestr(now);

 if SaveMe
     if ~exist('Outputs', 'dir'); mkdir('Outputs'); end
     DiaryStr = sprintf('DiarySingRunt%d.log',trial);
     diary(DiaryStr);
     runfile  = 'IsoRunLog.log';
     rlId = fopen(runfile,'a+');
%      fprintf(rlId, 'Isotropic Hard Rod \n');
    
 end

% Print what you are doing

fprintf('Isotropic Hard Rod \n')
fprintf('Start: %s\n', DateTimeStart)

%%%%%% Turn on/off interactions%%%%%%%%%
Interactions = 1;
SaveMe       = 1;
MakeMovies   = 0; % Movies won't run if save is zero
MakeOP       = 0;

%%%%%%%%%%%%% Box and Rod Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx      = 32;
Ny      = 32;
Nm      = 32;

%%%%%%%%% Initial density parameters%%%%%%%%%%%%%%%%%%
% Dimensionless  scaled concentration bc > 1.501 or bc < 1.499 if
% perturbing about equilbrum
bc      = 1.15;
L_rod   = 1;                  % Length of the rods
Lx      = 10*L_rod;               % Box length
Ly      = 10*L_rod;               % Box length
AspctRt = 8;                  % L / W
vD      = 0;                  %Driving velocity

%%%%%%%%%%%%%%%Time recording %%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t     = 1e-3; %time step
t_record    = 1e-2; %time interval for recording dynamics
t_tot       = 1e-2;   %total time
ss_epsilon  = 1e-8;                          %steady state condition

% The number of k-modes above and below k = 0 added as a perturbation
% Type of Inital Condition
% 0: Plane wave perturbation over an isotropic distribution
% 1: Plane wave perturbation over the equilibrium distribution
% 2: Plane wave perturbation over a nematic distribution
% 3: Load an equilbrium distribution
% 4: Seperate plane waves over an isotropic distribution (non-sensical)
% 5: Seperate plane waves over an nematic distribution (non-sensical)
% 6: A gaussian initial condition
IntCond     = 1;
NumModesX   = 4;
NumModesY   = 4;
NumModesM   = 4;
% Weight of the spatial sinusoidal perturbation. %
% Perturbations added to rho(i,j,k) = 1. Must be small
WeightPos   = 1e-3;
WeightAng   = 1e-3;
Random      = 0;       % Random perturbation coeffs

% Stepping method
% 0: AB1
% 1: AB2
% 2: HAB1
% 3: HAB2
% 4: BHAB1
% 5: BHAB2
% 6: Exponential Euler

StepMeth = 6; 

%%%%%%%%%%%%%%%%%%%%% Physical Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmp      = 1;            % Temperature
% mobility
Mob_same = 1;
Mob_pos  = Mob_same;
Mob_rot  = Mob_same;

[IntConcStr] =  IntDenNameWriter( IntCond );

if  strcmp( IntConcStr,'PlaneWaveEq' ) || strcmp( IntConcStr,'SepPWeq' )
    if 1.499 < bc && bc < 1.501
        bc = 1.502;
    end
end

if strcmp( IntConcStr,'Loaded' )
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
if vD  == 0; Drive = 0; else Drive = 1;end

% Parameter vectors
Paramtmp = [trial Interactions Drive StepMeth IntCond MakeOP MakeMovies SaveMe Nx Ny Nm...
    Lx Ly L_rod Tmp Norm WeightPos WeightAng Random ...
    NumModesX NumModesY NumModesM bc c Mob_pos Mob_rot vD];
Timetmp  = [delta_t t_record t_tot ss_epsilon];

% Make the output directory string and input file
FileDir = ...
    sprintf('HrIDC_N%i%i%i_bc%.2f_Int%i_v%.1f_IC%dt%ism%d',...
    Nx,Ny,Nm,bc,Interactions,vD,IntCond,trial,StepMeth);
FileInpt = ...
    sprintf('Inpt_N%i%i%i_bc%.2f_Int%i_v%.1f_IC%dt%i.txt', ...
    Nx,Ny,Nm,bc,Interactions,vD,IntCond,...
    trial);

Where2SavePath   = sprintf('%s/%s/%s',pwd,'Outputs',FileDir);

[fid,msg] = fopen(FileInpt,'w+');  % Note the 'wt' for writing in text mode

if fid == -1
    error( 'Failed to open %s: %s', FileInpt, msg);
end

fprintf(fid,'%s\n%s\n%s\n%s\n',FileDir,Where2SavePath,IntConcStr,IntDenName);
fprintf(fid,'Param\t');
fprintf(fid,'%e\t',Paramtmp);
fprintf(fid,'\n');
fprintf(fid,'Time\t');
fprintf(fid,'%e\t',Timetmp);

fclose(fid);

% keyboard
if Move == 1
    movefile(FileInpt,'/home/mws/Documents/Research/BG/DDFT/Inputs')
end

% keyboard
if Run == 1
    tic
    %     keyboard
    [DenFinal, DenFTFinal, GridObj, ParamObj,TimeObj,...
        DidIBreak,SteadyState,MaxReldRho] = ...
        HR2DrotMainIdC(FileInpt);
    toc
    disp('Params');disp(ParamObj);disp('Time');disp(TimeObj);      
    fprintf('Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
        DidIBreak,SteadyState,MaxReldRho)
    DateTimeEnd =  datestr(now);
    fprintf('End: %s\n\n', DateTimeEnd);
    
    if SaveMe
        diary off; 
        if ~exist(Where2SavePath, 'dir'); mkdir(Where2SavePath); end
        movefile('*.mat', Where2SavePath)
        movefile('*.txt', Where2SavePath)
        if MakeMovies
            movefile('*.avi', Where2SavePath)
        end
        if MakeOP
            movefile('*.fig', Where2SavePath)
            movefile('*.jpg', Where2SavePath)
        end         
        rlId = fopen(runfile,'a+');
        fprintf(rlId,'%s\n', FileInpt(1:end-4));
        fprintf(rlId,'Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
        DidIBreak,SteadyState,MaxReldRho);
        fprintf(rlId, 'Start: %s\n', DateTimeStart);
        fprintf(rlId,'End: %s\n\n', DateTimeEnd);
        fclose('all');
   end

end





