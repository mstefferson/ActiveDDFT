% Date
dateTime =  datestr(now);
fprintf('Starting RunHardRod: %s\n', dateTime);

% Add Subroutine path
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );

% Make Output Directories
if ~exist('runfiles', 'dir'); mkdir('./runfiles'); end;
if ~exist('runOPfiles', 'dir'); mkdir('./runOPfiles'); end;
if ~exist('analyzedfiles', 'dir'); mkdir('./analyzedfiles'); end;

% Grab initial parameters
if exist('Params.mat','file') == 0;
  if exist('initParams.m','file') == 0;
    cpParams();
  end;
  initParams();
end
load Params.mat;
movefile('Params.mat', 'ParamsRunning.mat');

% Copy the master parameter list to ParamObj
%ParamObj = ParamMaster;
systemObj = systemMaster;
particleObj = particleMaster;
rhoInit  = rhoInitMaster;
flags    = flagMaster;
runObj  = runMaster;
anisoDiffFlag = flags.AnisoDiff;
trial = runObj.trialID;

% Print what you are doing
if anisoDiffFlag  == 1;
  fprintf('Anisotropic Hard Rod \n')
else
  fprintf('Isotropic Hard Rod \n')
end

% Fix the time
fprintf('Making time obj\n');
[timeObj]= ...
  TimeStepRecMaker(timeMaster.dt,timeMaster.t_tot,...
    timeMaster.t_rec,timeMaster.t_write);
timeObj.ss_epsilon = timeMaster.ss_epsilon;
timeObj.amp_cutoff = timeMaster.amp_cutoff;
fprintf('Finished time obj\n');

% Display everythin
disp(flags); disp(particleObj); disp(systemObj); disp(timeObj); disp(rhoInit);

% Make paramMat
fprintf('Building parameter mat \n');
[paramMat, numRuns] = MakeParamMat( systemObj, particleObj, runObj, rhoInit, flags );
fprintf('Executing %d runs \n\n', numRuns);

paramvec = zeros(numRuns,1);
% For some reason, param_mat gets "sliced". Create vectors to get arround
paramNx  = paramMat(:,1); paramNy  = paramMat(:,2);
paramNm  = paramMat(:,3); paramLx  = paramMat(:,4);
paramLy  = paramMat(:,5); paramvD  = paramMat(:,6);
parambc  = paramMat(:,7); paramIC  = paramMat(:,8);
paramSM  = paramMat(:,9); paramrun = paramMat(:,10);

% Loops over all run
fprintf('Starting loop over runs\n');
if numRuns > 1
  parfor ii = 1:numRuns
    % Assign parameters
    paramvec = [ paramNx(ii) paramNy(ii) paramNm(ii) paramLx(ii) ...
      paramLy(ii) paramvD(ii) parambc(ii) paramIC(ii)...
      paramSM(ii) paramrun(ii)];
    
    % Name the file
    filename = [ 'Hr_Ani' num2str( anisoDiffFlag ) ...
      '_N' num2str( paramNx(ii) ) num2str( paramNy(ii) ) num2str( paramNm(ii) )  ...
      '_Lx' num2str( paramLx(ii) ) 'Ly' num2str( paramLy(ii) )...
      '_vD' num2str( paramvD(ii) ) '_bc' num2str( parambc(ii) ) ...
      '_IC' num2str( paramIC(ii) ) '_SM' num2str( paramSM(ii) ) ...
      '_t' num2str( trial ) '.' num2str( paramrun(ii) ) '.mat' ];
    
    fprintf('\nStarting %s \n', filename);
    [DenRecObj] = HR2DrotMain( filename, paramvec, systemObj, particleObj,...
       runObj, timeObj, rhoInit, flags );
    fprintf('Finished %s \n', filename);
  end
else
  paramvec = [ paramNx(1) paramNy(1) paramNm(1) paramLx(1) ...
    paramLy(1) paramvD(1) parambc(1) paramIC(1)...
    paramSM(1) paramrun(1)];
  
  % Name the file
  filename = [ 'Hr_Ani' num2str( anisoDiffFlag ) ...
    '_N' num2str( paramNx(1) ) num2str( paramNy(1) ) num2str( paramNm(1) )  ...
    '_Lx' num2str( paramLx(1) ) 'Ly' num2str( paramLy(1) )...
    '_vD' num2str( paramvD(1) ) '_bc' num2str( parambc(1) ) ...
    '_IC' num2str( paramIC(1) ) '_SM' num2str( paramSM(1) ) ...
    '_t' num2str( trial ) '.' num2str( paramrun(1) ) '.mat' ];
  
  fprintf('\nStarting %s \n', filename);
  [DenRecObj] = HR2DrotMain( filename, paramvec, systemObj, particleObj,...
    runObj, timeObj, rhoInit, flags );
  fprintf('Finished %s \n', filename);
end


dateTime =  datestr(now);
movefile('ParamsRunning.mat', 'ParamsFinished.mat');
fprintf('Finished RunHardRod: %s\n', dateTime);

