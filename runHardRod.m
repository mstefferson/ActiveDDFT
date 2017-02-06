% Executeable for HardRod
% Date
function denRecObj = runHardRod()
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
if exist('Params.mat','file') == 0
  if exist('initParams.m','file') == 0
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
diagOp = flags.DiagLop;
trial = runObj.trialID;
% grab particle types and interactions
ptype = ['_' particleObj.type];
if isempty(particleObj.interHb); interHb = '';
else; interHb = ['_' particleObj.interHb]; end
if isempty(particleObj.interLr); interLr = '';
else; interLr = ['_' particleObj.interLr]; end
if isempty(particleObj.externalPot); externalPot = '';
else; externalPot = ['_' particleObj.externalPot]; end
% Print what you are doing
if diagOp  == 1
  fprintf('Diagonal operator (cube) \n')
else
  fprintf('Off-diagon operator (expokit) \n')
end
% Scramble seed if you want
if flags.rndStrtUpSeed
  rng('shuffle');
  fprintf('Shuffling startup seed\n')
else
  fprintf('Using MATLABs original seed\n');
end
% Fix the time
fprintf('Making time obj\n');
[timeObj]= ...
  TimeStepRecMaker(timeMaster.dt,timeMaster.t_tot,...
  timeMaster.t_rec,timeMaster.t_write);
timeObj.ss_epsilon = timeMaster.ss_epsilon;
timeObj.ss_epsilon_dt = timeMaster.ss_epsilon_dt;
fprintf('Finished time obj\n');
% Display everythin
disp(runObj); disp(flags); disp(particleObj); disp(systemObj); disp(timeObj); disp(rhoInit);
% Make paramMat
fprintf('Building parameter mat \n');
[paramMat, numRuns] = MakeParamMat( systemObj, particleObj, runObj, rhoInit, flags );
fprintf('Executing %d runs \n\n', numRuns);
% For some reason, param_mat gets "sliced". Create vectors to get arround
paramn1  = paramMat(:,1); paramn2  = paramMat(:,2);
paramn3  = paramMat(:,3); paraml1  = paramMat(:,4);
paraml2  = paramMat(:,5); paramvD  = paramMat(:,6);
parambc  = paramMat(:,7); paramIC  = paramMat(:,8);
paramSM  = paramMat(:,9); paramrun = paramMat(:,10);
% Loops over all run
fprintf('Starting loop over runs\n');
ticID = tic;
if numRuns > 1
  parobj = gcp;
  fprintf('I have hired %d workers\n',parobj.NumWorkers);
  parfor ii = 1:numRuns
    % Assign parameters
    paramvec = [ paramn1(ii) paramn2(ii) paramn3(ii) paraml1(ii) ...
      paraml2(ii) paramvD(ii) parambc(ii) paramIC(ii)...
      paramSM(ii) paramrun(ii)];
    % Name the file
    filename = [ 'Hr' ptype, interHb, interLr,  externalPot, ...
      '_diag' num2str( diagOp ) ...
      '_N' num2str( paramn1(ii) ) num2str( paramn2(ii) ) num2str( paramn3(ii) )  ...
      '_ls' num2str( paraml1(ii) ) num2str( paraml2(ii) )...
      '_bc' num2str( parambc(ii) ) '_vD' num2str( paramvD(ii) ) ...
      '_IC' num2str( paramIC(ii) ) '_SM' num2str( paramSM(ii) ) ...
      '_t' num2str( trial,'%.2d' ) '.' num2str( paramrun(ii), '%.2d' ) '.mat' ];
    fprintf('\nStarting %s \n', filename);
    [denRecObj] = ddftMain( filename, paramvec, systemObj, particleObj,...
      runObj, timeObj, rhoInit, flags );
    fprintf('Finished %s \n', filename);
  end
else
  paramvec = [ paramn1(1) paramn2(1) paramn3(1) paraml1(1) ...
    paraml2(1) paramvD(1) parambc(1) paramIC(1)...
    paramSM(1) paramrun(1)];
  % Name the file
  filename = [ 'Hr' ptype, interHb, interLr,  externalPot, ...
    '_diag' num2str( diagOp ) ...
    '_N' num2str( paramn1(1) ) num2str( paramn2(1) ) num2str( paramn3(1) )  ...
    '_ls' num2str( paraml1(1) )  num2str( paraml2(1) )...
    '_bc' num2str( parambc(1) ) '_vD' num2str( paramvD(1) ) ...
    '_IC' num2str( paramIC(1) ) '_SM' num2str( paramSM(1) ) ...
    '_t' num2str( trial,'%.2d' ) '.' num2str( paramrun(1), '%.2d' ) '.mat' ];
  fprintf('\nStarting %s \n', filename);
  [denRecObj] = ddftMain( filename, paramvec, systemObj, particleObj,...
    runObj, timeObj, rhoInit, flags );
  fprintf('Finished %s \n', filename);
end
runTime = toc(ticID);
dateTime =  datestr(now);
movefile('ParamsRunning.mat', 'ParamsFinished.mat');
runHr = floor( runTime / 3600); runTime = runTime - runHr*3600;
runMin = floor( runTime / 60);  runTime = runTime - runMin*60;
runSec = floor(runTime);
fprintf('RunTime: %.2d:%.2d:%.2d (hr:min:sec)\n', runHr, runMin,runSec);
fprintf('Finished RunHardRod: %s\n', dateTime);
% check for log file
logFile = dir('*.out');
if ~isempty(logFile)
  if length(logFile) == 1
    fprintf('One log file found. Copying it to save dir.\n')
  else
    fprintf('Too many files. Copying all to save dir.\n')
  end
  for ii = 1:length(logFile)
    newLog = [ filename(1:end-4) '_l' num2str(ii,'%.2d') '.out' ];
    movefile(logFile(ii).name, newLog );
    movefile( newLog, denRecObj.dirName );
  end
else
  fprintf('No log file found\n')
end