
% Date
DateTimeStart =  datestr(now);
fprintf('Starting RunHardRod: %s\n', DateTimeStart);

% Add Subroutine path
CurrentDir = pwd;
addpath( genpath( [CurrentDir '/Subroutines'] ) );

% Make Output Directories
if ~exist('runfiles', 'dir'); mkdir('./runfiles'); end;
if ~exist('runOPfiles', 'dir'); mkdir('./runOPfiles'); end;
if ~exist('analyzedfiles', 'dir'); mkdir('./analyzedfiles'); end;

% Grab initial parameters
if exist('Params.mat','file') == 0;
  if exist('InitParams.m','file') == 0;
    cpParams
  end;
  InitParams
end
load Params.mat;

% Copy the master parameter list to ParamObj
ParamObj = ParamMaster;
RhoInit  = RhoInitMaster;
Flags    = FlagMaster;

% Print what you are doing
if Flags.AnisoDiff  == 1;
  fprintf('Anisotropic Hard Rod \n')
else
  fprintf('Isotropic Hard Rod \n')
end

% Fix the time
[TimeObj.t_tot,TimeObj.N_time,TimeObj.t_rec,TimeObj.N_rec,TimeObj.N_count]= ...
        TimeStepRecMaker(TimeObj.delta_t,TimeObj.t_tot,TimeObj.t_record);

% Make paramMat 
[paramMat, numRuns] = MakeParamMat( ParamObj, RhoInit, Flags );

% For some reason, param_mat gets "sliced". Create vectors to get arround
paramNx  = paramMat(:,1); paramNy  = paramMat(:,2); paramNm  = paramMat(:,3);
paramLx  = paramMat(:,4); paramLy  = paramMat(:,5); paramvD  = paramMat(:,6);
parambc  = paramMat(:,7); paramIC  = paramMat(:,8); paramSM  = paramMat(:,9);
paramrun = paramMat(:,10);

% Loops over all run

for ii = 1:numRuns
  % Assign parameters
  ParamObj.Nx = paramNx(ii);
  ParamObj.Ny = paramNy(ii);
  ParamObj.Nm = paramNm(ii);
  ParamObj.Lx = paramLx(ii);
  ParamObj.Ly = paramLy(ii);
  ParamObj.vD = paramvD(ii);
  ParamObj.bc = parambc(ii);
  RhoInit.IntCond = paramIC(ii);
  Flags.StepMeth = paramSM(ii);
  ParamObj.runID = paramrun(ii);
  ParamObj.Norm  = ParamObj.bc / ( ParamObj.L_rod ^ 2 / pi ) *... % number of particles
    ParamObj.Lx * ParamObj.Ly;
    ParamObj.c = ParamObj.bc ./ ParamObj.b;
  if Flags.SquareBox
    ParamObj.L_box = ParamObj.Lx;
  else
    ParamObj.L_box = [ParamObj.Lx ParamObj.Ly];
  end
  
  % Name the file
  filename = [ 'Hr_Ani' num2str( Flags.AnisoDiff ) ...
    '_N' num2str( ParamObj.Nx ) num2str( ParamObj.Ny ) num2str( ParamObj.Nm )  ...
    '_Lx' num2str( ParamObj.Lx ) 'Ly' num2str( ParamObj.Ly )...
    '_vD' num2str( ParamObj.vD ) '_bc' num2str( ParamObj.bc ) ...
    '_IC' num2str( RhoInit.IntCond ) '_SM' num2str( Flags.StepMeth ) ...
    '_t' num2str( ParamObj.trial ) '.' num2str( ParamObj.runID ) '.mat' ];
  
  disp(filename);
  disp(Flags);
  disp(ParamObj);
  disp(RhoInit);
  disp(TimeObj);

  [DidIBreak,SteadyState,MaxReldRho] = ...
      HR2DrotMain( filename, ParamObj, TimeObj, RhoInit, Flags );
      
  
end

delete Params.mat
%%
%if  strcmp( IntConcStr,'PlaneWaveEq' ) || strcmp( IntConcStr,'SepPWeq' )
%if 1.499 < bc && bc < 1.501
%bc = 1.502;
%end
%end


%% Make the output directory string and input file
%if SaveMe
%if ~exist('Outputs', 'dir'); mkdir('Outputs'); end;
%DiaryStr = sprintf('DiarySingRunt%d.log',trial);
%diary(DiaryStr);
%runfile  = 'AnisRunLog.log';
%rlId = fopen(runfile,'a+');
%%      fprintf(rlId, 'Anisotropic Hard Rod Hard Rod \n');
%end

%FileDir = ...
%sprintf('HrAnD%i%i%i_bc%.2f_Int%i_v%.1f_IC%dsm%dt%d',...
%Nx,Ny,Nm,bc,Interactions,vD,IntCond,StepMeth,trial);
%FileInpt = ...
%sprintf('InptAnD_N%i%i%i_bc%.2f_Int%i_v%.1f_IC%dsm%dt%d.txt', ...
%Nx,Ny,Nm,bc,Interactions,vD,IntCond,StepMeth,trial);

%Where2SavePath    = sprintf('%s/%s/%s',pwd,'Outputs',FileDir);

%% Figure out Total Runs
%tic
%%     keyboard
%[DenFinal, DenFTFinal, GridObj, ParamMaster,TimeObj,...
%DidIBreak,SteadyState,MaxReldRho] = ...
%HR2DrotMain(FileInpt);
%toc
%disp('Params');disp(ParamMaster);disp('Time');disp(TimeObj);
%fprintf('Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
%DidIBreak,SteadyState,MaxReldRho)
%DateTimeEnd =  datestr(now);
%fprintf('End: %s\n\n', DateTimeEnd);

%if SaveMe
%diary off;
%if ~exist(Where2SavePath, 'dir'); mkdir(Where2SavePath); end
%movefile('*.mat', Where2SavePath)
%movefile('*.txt', Where2SavePath)
%movefile('Diary*', Where2SavePath)
%if MakeMovies
%movefile('*.avi', Where2SavePath)
%end
%if MakeOP
%movefile('*.fig', Where2SavePath)
%movefile('*.jpg', Where2SavePath)
%end

%rlId = fopen(runfile,'a+');
%fprintf(rlId,'%s\n', FileInpt(1:end-4));
%fprintf(rlId,'( dt t_rec t_tot ss_eps ): ( %.2e %.2e %.2e %.2e )\n',...
%Timetmp);
%fprintf(rlId,'( Nx, Ny, Nm ): ( %d, %d, %d )\n', Nx, Ny, Nm);
%fprintf(rlId,'( Lx, Ly, Lrod ) = ( %.1f, %.1f, %.1f )\n',...
%Lx, Ly, L_rod);
%fprintf(rlId,...
%'( bc, vD, Mob_par, Mob_perp, Mob_rot ) = ( %.2f, %.2f, %.1f, %.1f, %.1f )\n',...
%bc, vD, Mob_par, Mob_perp, Mob_rot);
%fprintf(rlId,'( MdX, MdY, MdM, rand) = ( %d, %d, %d, %d )\n',...
%NumModesX, NumModesY, NumModesM, RandomAmp);
%fprintf(rlId,'Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
%DidIBreak,SteadyState,MaxReldRho);
%fprintf(rlId, 'Start: %s\n', DateTimeStart);
%fprintf(rlId,'End: %s\n\n', DateTimeEnd);
%fclose('all');
%else
%remove('*.txt')
%end

%end

