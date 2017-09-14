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
% Copy the master parameter list to ParamObj
systemObj = systemMaster;
particleObj = particleMaster;
timeObj = timeMaster;
rhoInit  = rhoInitMaster;
flags    = flagMaster;
runObj  = runMaster;
diagOp = flags.DiagLop;
trial = runObj.trialID;
% Fix things
% Change odd gridspacings to even unless it's one.
if systemObj.n1 == 1
  systemObj.l1 = 1;
else
  systemObj.n1 = systemObj.n1 + mod( systemObj.n1, 2 );
end
if systemObj.n2 == 1
  systemObj.l2 = 1;
else
  systemObj.n2 = systemObj.n2 + mod( systemObj.n2, 2 );
end
if systemObj.n3 == 1
  systemObj.l3 = 1;
else
  systemObj.n3 = systemObj.n3 + mod( systemObj.n3, 2 );
end
% Fix Ls if we want the box to be square
if flags.SquareBox == 1
  systemObj.L_box = unique( [systemObj.l1 systemObj.l2] );
  systemObj.l1 = systemObj.L_box;
  systemObj.l2 = systemObj.L_box;
end
% Fix l1 is we want all Ns to be the same
if flags.AllNsSame == 1
  if systemObj.n3 == 1
    Nvec = unique( [systemObj.n1 systemObj.n2] );
    systemObj.n1 = Nvec;  systemObj.n2 = Nvec;
  else
    Nvec = unique( [systemObj.n1 systemObj.n2 systemObj.n3] );
    systemObj.n1 = Nvec;  systemObj.n2 = Nvec;   systemObj.n3 = Nvec;
  end
end
% Make OP if making movies
if flags.MakeMovies == 1; flags.MakeOP = 1; end % if make movie, make OP first
%Currently, you must save
if flags.MakeOP && flags.SaveMe == 0
  fprintf('Turning on saving, you must be saving to make OPs (due to matfile)\n');
  flags.SaveMe = 1;
end
if particleObj.vD  == 0; flags.Drive = 0; else flags.Drive = 1;end
% fix rhoInit
rhoInitObj = checkRhoInit( rhoInit, systemObj );
% Get particle mobility
[particleObj, systemObj] = ...
  particleInit( particleObj, systemObj, flags.DiagLop);
% grab particle types and interactions
ptype = ['_' particleObj.type];
if isempty(particleObj.interHb); interHb = '';
else; interHb = ['_' particleObj.interHb]; end
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
% external potentials
fprintf('Building external potential \n' )
nVec = [systemObj.n1 systemObj.n2 systemObj.n3];
[ particleObj.externalV] = ...
  checkExtVDim( particleObj.externalV, nVec );
[ particleObj.interactLrV] = ...
  checkInteractVDim( particleObj.interactLrV, systemObj.n3 );
[externVObj] = potRunManager( particleObj.externalV, 0 );
% mean field interactions
[interactLrVObj] = potRunManager( particleObj.interactLrV, 1 );
% Make paramMat
fprintf('Building parameter mat \n');
[paramMat, numRuns] = ...
  MakeParamMat( systemObj, particleObj, runObj, externVObj.inds, ...
  interactLrVObj.inds, flags );
fprintf('Executing %d runs \n\n', numRuns);
% Display everythin
disp(runObj); disp(flags); disp(particleObj); disp(systemObj); disp(timeObj); disp(rhoInitObj);
% For some reason, param_mat gets "sliced". Create vectors to get arround
paramn1      = paramMat(1,:);  paramn2 = paramMat(2,:);
paramn3      = paramMat(3,:);  paraml1 = paramMat(4,:);
paraml2      = paramMat(5,:);  paramvD = paramMat(6,:);
parambc      = paramMat(7,:);  paramSM = paramMat(8,:); 
paramRun     = paramMat(9,:); paramInteractInds = paramMat(10,:);
paramExtInds = paramMat(11,:);
% potentials cells
extParam = externVObj.param;
extStr = externVObj.str;
interactParam = interactLrVObj.param;
interactStr = interactLrVObj.str;
% rhoInit str
initStr = rhoInit.intCond{1};
keyboard
% Loops over all run
fprintf('Starting loop over runs\n');
ticID = tic;
if numRuns > 1
  %parobj = gcp;
  %fprintf('I have hired %d workers\n',parobj.NumWorkers);
  for ii = 1:numRuns
    % Assign parameters
    paramvec = [ paramn1(ii) paramn2(ii) paramn3(ii) paraml1(ii) ...
      paraml2(ii) paramvD(ii) parambc(ii) paramIC(ii)...
      paramSM(ii)  paramRun(ii)];
    % Name the file
    filename = [ 'Hr' ptype, interHb, ...
    interactStr{paramInteractInds(ii)},  extStr{paramExtInds(ii)}, ...
      '_diag' num2str( diagOp ) ...
      '_N' num2str( paramn1(ii) ) num2str( paramn2(ii) ) num2str( paramn3(ii) )  ...
      '_ls' num2str( paraml1(ii) ) num2str( paraml2(ii) )...
      '_bc' num2str( parambc(ii), '%.2f' ) '_vD' num2str( paramvD(ii), '%.3g' ) ...
      '_IC' num2str( paramIC(ii), '%d' ) '_SM' num2str( paramSM(ii), '%d'  ) ...
      '_t' num2str( trial,'%.2d' ) '.' num2str( paramRun(ii), '%.2d' ) '.mat' ];
    fprintf('\nStarting %s \n', filename);
    [denRecObj] = ddftMain( filename, paramvec, systemObj, particleObj,...
      runObj, timeObj, rhoInitObj, flags, ...
      interactParam{paramInteractInds(ii)},extParam{paramExtInds(ii)} );
    fprintf('Finished %s \n', filename);
  end
else
  ii = 1;
    % Assign parameters
    paramvec = [ paramn1(ii) paramn2(ii) paramn3(ii) paraml1(ii) ...
      paraml2(ii) paramvD(ii) parambc(ii) paramIC(ii)...
      paramSM(ii)  paramRun(ii)];
    % Name the file
    filename = [ 'Hr' ptype, interHb, ...
    interactStr{paramInteractInds(ii)},  extStr{paramExtInds(ii)}, ...
      '_diag' num2str( diagOp ) ...
      '_N' num2str( paramn1(ii) ) num2str( paramn2(ii) ) num2str( paramn3(ii) )  ...
      '_ls' num2str( paraml1(ii) ) num2str( paraml2(ii) )...
      '_bc' num2str( parambc(ii), '%.2f' ) '_vD' num2str( paramvD(ii), '%.3g' ) ...
      '_IC' num2str( paramIC(ii), '%d' ) '_SM' num2str( paramSM(ii), '%d'  ) ...
      '_t' num2str( trial,'%.2d' ) '.' num2str( paramRun(ii), '%.2d' ) '.mat' ];
    fprintf('\nStarting %s \n', filename);
    [denRecObj] = ddftMain( filename, paramvec, systemObj, particleObj,...
      runObj, timeObj, rhoInitObj, flags, ...
      interactParam{paramInteractInds(ii)},extParam{paramExtInds(ii)} );
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
