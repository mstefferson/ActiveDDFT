% Executeable for HardRod
% Date
function denRecObj = runHardRod()
dateTime =  datestr(now);
fprintf('Starting RunHardRod: %s\n', dateTime);
% Add Subroutine path
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% Make Output Directories
if ~exist('runfiles', 'dir'); mkdir('./runfiles'); end
if ~exist('runOPfiles', 'dir'); mkdir('./runOPfiles'); end
if ~exist('analyzedfiles', 'dir'); mkdir('./analyzedfiles'); end
% Grab initial parameters
if exist('Params.mat','file') == 0
  if exist('initParams.m','file') == 0
    cpParams();
  end
  fprintf('No Params.mat found. Running initParams\n');
  initParams();
end
load Params.mat;
fprintf('Params locked and loaded\n');
movefile('Params.mat', 'ParamsRunning.mat');
% Fix things
% Change odd gridspacings to even unless it's one.
if systemMaster.n1 == 1
  systemMaster.l1 = 1;
  rhoInitMaster.NumModesX = 0;
else
  systemMaster.n1 = systemMaster.n1 + mod( systemMaster.n1, 2 );
end
if systemMaster.n2 == 1
  systemMaster.l2 = 1;
  rhoInitMaster.NumModesY = 0;
else
  systemMaster.n2 = systemMaster.n2 + mod( systemMaster.n2, 2 );
end
if systemMaster.n3 == 1
  systemMaster.l3 = 1;
  rhoInitMaster.NumModesM = 0;
else
  systemMaster.n3 = systemMaster.n3 + mod( systemMaster.n3, 2 );
end
% Don't perturb more more than you are allowed to
if( rhoInitMaster.NumModesX >= systemMaster.n1 / 2 )
  rhoInitMaster.NumModesX = floor(systemMaster.n1 / 2) - 1;
end
if( rhoInitMaster.NumModesY >= systemMaster.n2 / 2 )
  rhoInitMaster.NumModesY = floor(systemMaster.n2 / 2) - 1;
end
if( rhoInitMaster.NumModesM >= systemMaster.n3 / 2 )
  rhoInitMaster.NumModesM = floor(systemMaster.n3 / 2) - 1;
end
%  Make sure variance isn't zero if doing polar
if rhoInitMaster.gP(2) == 0;  rhoInitMaster.gP(2) = min(systemMaster.l1)/2; end
if rhoInitMaster.gP(5) == 0;  rhoInitMaster.gP(2) = min(systemMaster.l1)/2; end
if rhoInitMaster.gP(8) == 0; rhoInitMaster.gP(8) = min(systemMaster.l1)/2; end
% Scale ss_epsilon by delta_t. Equivalent to checking d rho /dt has reached
% steady state instead of d rho
timeMaster.ss_epsilon_dt = timeMaster.ss_epsilon .* timeMaster.dt;
% Fix Ls if we want the box to be square
if flagMaster.SquareBox == 1
  systemMaster.L_box = unique( [systemMaster.l1 systemMaster.l2] );
  systemMaster.l1 = systemMaster.L_box;
  systemMaster.l2 = systemMaster.L_box;
end
% Fix l1 is we want all Ns to be the same
if flagMaster.AllNsSame == 1
  if systemMaster.n3 == 1
    Nvec = unique( [systemMaster.n1 systemMaster.n2] );
    systemMaster.n1 = Nvec;  systemMaster.n2 = Nvec;
  else
  Nvec = unique( [systemMaster.n1 systemMaster.n2 systemMaster.n3] );
  systemMaster.n1 = Nvec;  systemMaster.n2 = Nvec;   systemMaster.n3 = Nvec;
  end
end
% Make OP if making movies
if flagMaster.MakeMovies == 1; flagMaster.MakeOP = 1; end % if make movie, make OP first
%Currently, you must save
if flagMaster.MakeOP && flagMaster.SaveMe == 0
  fprintf('Turning on saving, you must be saving to make OPs (due to matfile)\n');
  flagMaster.SaveMe = 1;
end
if particleMaster.vD  == 0; flagMaster.Drive = 0; else flagMaster.Drive = 1;end
% Get particle mobility
[particleMaster, systemMaster] = ...
  particleInit( particleMaster, systemMaster, flagMaster.DiagLop);
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
paramn1  = paramMat(1,:); paramn2  = paramMat(2,:);
paramn3  = paramMat(3,:); paraml1  = paramMat(4,:);
paraml2  = paramMat(5,:); paramvD  = paramMat(6,:);
parambc  = paramMat(7,:); paramIC  = paramMat(8,:);
paramSM  = paramMat(9,:); paramrun = paramMat(10,:);
paramLrL1  = paramMat(11,:); paramLrL2 = paramMat(12,:);
paramLrE1  = paramMat(13,:); paramLrE2 = paramMat(14,:);
% Get interaction strings
lrStr = interLr;
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
      paramSM(ii) paramrun(ii) paramLrL1(ii) paramLrL2(ii) paramLrE1(ii) paramLrE2(ii) ];
    % long range string
    if strcmp( interLr, '_softshoulder' )
      lrStrTemp = [ interLr '_1il' num2str( paramLrL1(ii), '%.2f' ) '_2il' ...
        num2str( paramLrL2(ii), '%.2f' )...
        '_1ie' num2str( paramLrE1(ii), '%.2f' )  '_2ie' ...
        num2str( paramLrE2(ii), '%.2f' ) ];
    else
      lrStrTemp = lrStr;
    end
    % Name the file
    filename = [ 'Hr' ptype, interHb, lrStrTemp,  externalPot, ...
      '_diag' num2str( diagOp ) ...
      '_N' num2str( paramn1(ii) ) num2str( paramn2(ii) ) num2str( paramn3(ii) )  ...
      '_ls' num2str( paraml1(ii) ) num2str( paraml2(ii) )...
      '_bc' num2str( parambc(ii), '%.2f' ) '_vD' num2str( paramvD(ii), '%.3g' ) ...
      '_IC' num2str( paramIC(ii), '%d' ) '_SM' num2str( paramSM(ii), '%d'  ) ...
      '_t' num2str( trial,'%.2d' ) '.' num2str( paramrun(ii), '%.2d' ) '.mat' ];
    fprintf('\nStarting %s \n', filename);
   [denRecObj] = ddftMain( filename, paramvec, systemObj, particleObj,...
      runObj, timeObj, rhoInit, flags );
    fprintf('Finished %s \n', filename);
  end
else
  ii = 1;
  % Assign parameters
  paramvec = [ paramn1(ii) paramn2(ii) paramn3(ii) paraml1(ii) ...
    paraml2(ii) paramvD(ii) parambc(ii) paramIC(ii)...
    paramSM(ii) paramrun(ii) paramLrL1(ii) paramLrL2(ii) paramLrE1(ii) paramLrE2(ii) ];
    % long range string
    if strcmp( interLr, '_softshoulder' )
      lrStrTemp = [ interLr '_1il' num2str( paramLrL1(ii), '%.2f' ) '_2il' ...
        num2str( paramLrL2(ii), '%.2f' )...
        '_1ie' num2str( paramLrE1(ii), '%.2f' )  '_2ie' ...
        num2str( paramLrE2(ii), '%.2f' ) ];
    else
      lrStrTemp = lrStr;
    end
    % Name the file
    filename = [ 'Hr' ptype, interHb, lrStrTemp,  externalPot, ...
      '_diag' num2str( diagOp ) ...
      '_N' num2str( paramn1(ii) ) num2str( paramn2(ii) ) num2str( paramn3(ii) )  ...
      '_ls' num2str( paraml1(ii) ) num2str( paraml2(ii) )...
      '_bc' num2str( parambc(ii), '%.2f' ) '_vD' num2str( paramvD(ii), '%.3g' ) ...
      '_IC' num2str( paramIC(ii), '%d' ) '_SM' num2str( paramSM(ii), '%d'  ) ...
      '_t' num2str( trial,'%.2d' ) '.' num2str( paramrun(ii), '%.2d' ) '.mat' ];
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
