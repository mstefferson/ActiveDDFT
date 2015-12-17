% Input creater for HR2DrotMainDrIDCube
%
% All subroutines must be located in current directory
% Add path

CurrentDir = pwd;
addpath( genpath( CurrentDir) );

% Now can change number of grid points in the x, y, phi direction
Run  = 1; % Run main from here
Move = 0; % Move files to a nice location

%%%%%%%% Trial %%%%%%%%%%%%
trial    = 4;

%%%%%% Turn on/off interactions%%%%%%%%%
Interactions = 1;
SaveMe       = 1;
MakeMovies   = 1; % Movies won't run if save is zero
MakeOP       = 1;
%%%%%%%%%%%%% Box and Rod Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx      = 64;
Ny      = 64;
Nm      = 64;

L_rod   = 1;                  % Length of the rods
Lx      = 10*L_rod;               % Box length
Ly      = 10*L_rod;               % Box length
AspctRt = 8;                  % L / W
vD      = 0;                  %Driving velocity

%%%%%%%%%%%%%%%Time recording %%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t     = 1e-3; %time step
t_record    = 1; %time interval for recording dynamics
t_tot       = 50;   %total time
ss_epsilon  = 1e-8;                          %steady state condition

% The number of k-modes above and below k = 0 added as a perturbation
NumModesX   = 4;
NumModesY   = 4;
NumModesM   = 4;
WeightPos   = 1e-3;
WeightAng   = 1e-3;
Random      = 0;       % Random perturbation coeffs

%%%%%%%%%%%%%%%%%%%%% Physical Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmp      = 1;            % Temperature
% mobility
Mob_same = 1;
Mob_pos  = Mob_same;
Mob_rot  = Mob_same;

%%%%%%%%% Initial density parameters%%%%%%%%%%%%%%%%%%
% Weight of the spatial sinusoidal perturbation. %
% Perturbations added to rho(i,j,k) = 1. Must be small
% Dimensionless  scaled concentration bc > 1.501 or bc < 1.499 if
% perturbing about equilbrum
bc             = 1.0;


% Type of initial Condition
IntGauss   = 0;
IntPw      = 1;
IntSepPw   = 0;
IntEqPw    = 0;    % Distribution from perturbing equil. dist.
IntEqSepPw = 0;    % This is the wrong way to perturb.
IntLoad    = 0;
IntNemPw   = 0;    % Distribution from perturbing equil. dist.

% Save a string saying what you want
[IntDenType, IntDenIndicator] = ...
    IntDenIndicatorMaker(IntGauss,IntPw,IntSepPw,IntEqPw,...
    IntEqSepPw,IntLoad,IntNemPw);

if IntEqPw == 1 || IntSepPw == 1
    if 1.499 < bc && bc < 1.501
        bc = 1.502;
    end
end

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
if vD  == 0; Drive = 0; else Drive = 1;end

% Parameter vectors
Paramtmp = [trial Interactions Drive MakeOP MakeMovies SaveMe Nx Ny Nm...
    Lx Ly L_rod Tmp Norm WeightPos WeightAng Random ...
    NumModesX NumModesY NumModesM bc c Mob_pos Mob_rot vD];
Timetmp  = [delta_t t_record t_tot ss_epsilon];

% Make the output directory string and input file
FileDir = ...
    sprintf('HRdiffIDC_N%i%i%i_bc%.1f_Int%i_v%.1f_%st%i',...
    Nx,Ny,Nm,bc,Interactions,vD,IntDenIndicator,trial);
FileInpt = ...
    sprintf('Inpt_N%i%i%i_bc%.1f_Int%i_v%.1f_%st%i.txt', ...
    Nx,Ny,Nm,bc,Interactions,vD,IntDenIndicator,...
    trial);

Where2SavePath    = sprintf('%s/%s/%s',pwd,'Outputs',FileDir);

[fid,msg] = fopen(FileInpt,'w+');  % Note the 'wt' for writing in text mode

if fid == -1
    error( 'Failed to open %s: %s', FileInpt, msg);
end

fprintf(fid,'%s\n%s\n%s\n',FileDir,Where2SavePath,IntDenType);
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
        HR2DrotMainDrIDCube(FileInpt);

    if SaveMe
        mkdir Outputs
        DiaryStr = sprintf('DiarySingRunt%d.txt',trial);
        diary(DiaryStr);
        disp('Params');disp(ParamObj);disp('Time');disp(TimeObj);
        mkdir(Where2SavePath)
        movefile('*.mat', Where2SavePath)
        movefile('*.txt', Where2SavePath)
        if MakeMovies
            movefile('*.avi', Where2SavePath)
        end
        if MakeOP
            movefile('*.fig', Where2SavePath)
            movefile('*.jpg', Where2SavePath)
        end  
        
    end
    toc
    fprintf('Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
        DidIBreak,SteadyState,MaxReldRho)
    diary off
    %     cd /home/mws/Documents/MATLAB/Research/BG/DDFT/HRddft/Drive/IsoDiffCube
end





