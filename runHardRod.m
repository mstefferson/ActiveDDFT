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
flags.movie = movieFlagMaster;
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
if flags.movie.analysis == 1; flags.MakeOP = 1; end % if make movie, make OP first
if flags.MakeOP == 1; flags.SaveMe = 1; end
%Currently, you must save
if flags.MakeOP && flags.SaveMe == 0
  fprintf('Turning on saving, you must be saving to make OPs (due to matfile)\n');
  flags.SaveMe = 1;
end
if particleObj.vD  == 0; flags.Drive = 0; else flags.Drive = 1;end
% Get particle mobility
[particleObj, systemObj] = ...
  particleInit( particleObj, systemObj, flags.DiagLop);
% fix rhoInit
rhoInitObj = rhoInitManager( rhoInit, systemObj, particleObj );
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
particleObj.externalV = ...
  checkExtVDim( particleObj.externalV, nVec );
particleObj.interactLrV = ...
  checkInteractVDim( particleObj.interactLrV, systemObj.n3 );
externVObj = potRunManager( particleObj.externalV, 0 );
noiseObj = noiseRunManager( systemObj.noise, systemObj.n3 );
% mean field interactions
interactLrVObj = potRunManager( particleObj.interactLrV, 1 );
% Make paramMat
fprintf('Building parameter mat \n');
[paramMat, numRuns] = ...
  MakeParamMat( systemObj, particleObj, runObj, externVObj.inds, ...
  interactLrVObj.inds, noiseObj.inds, flags );
fprintf('Executing %d runs \n\n', numRuns);
% Display everythin
disp(runObj); disp(flags); disp(particleObj); disp(systemObj); disp(timeObj); disp(rhoInitObj);
% For some reason, param_mat gets "sliced". Create vectors to get arround
paramn1  = paramMat(1,:); paramn2  = paramMat(2,:);
paramn3  = paramMat(3,:); paraml1  = paramMat(4,:);
paraml2  = paramMat(5,:); paramvD  = paramMat(6,:);
parambc  = paramMat(7,:); paramSM  = paramMat(8,:); 
paramRun = paramMat(9,:);
paramInteractInds = paramMat(10,:);
paramExtInds = paramMat(11,:);
paramNoiseInds = paramMat(12,:);
% potentials cells and noise
extParam = externVObj.param;
extStr = externVObj.str;
noiseParam = noiseObj.param;
interactParam = interactLrVObj.param;
% set-up strs
interactStr = cell(1,numRuns);
extStr = cell(1,numRuns);
noiseStr = cell(1,numRuns);
% nlDr
nlDiff = [ particleObj.nlDiff{1} '_' ...
  num2str( particleObj.nlDiff{2} )...
  '_' num2str( particleObj.nlDiff{3}, '%d' ) ];
for ii = 1:numRuns
  interactStr{ii} = interactLrVObj.str{ paramInteractInds(ii) };
  extStr{ii} = externVObj.str{ paramExtInds(ii) };
  noiseStr{ii} = noiseObj.str{ paramNoiseInds(ii) };
end
% rhoInit str
initStr = rhoInitObj.fileStr;
% set-up parfor
if numRuns > 1 && flags.parforFlag
  parobj = gcp;
  numWorkers = parobj.NumWorkers;
  flags.movie.analysis = 0;
else
  numWorkers = 0;
end
% Loops over all run
denRecObj = cell(numRuns,1);
fprintf('Starting loop over runs\n');
ticID = tic;
parfor (ii = 1:numRuns, numWorkers)
  % Assign parameters
  paramvec = [ paramn1(ii) paramn2(ii) paramn3(ii) paraml1(ii) ...
    paraml2(ii) paramvD(ii) parambc(ii) paramSM(ii)  paramRun(ii)...
    paramInteractInds(ii) paramExtInds(ii) paramNoiseInds(ii) ];
  % Name the file
  filename = [ 'Hr' ptype, interHb, ...
    interactStr{ ii },  extStr{ ii }, ...
    '_diag' num2str( diagOp ) ...
    '_N' num2str( paramn1(ii) ) num2str( paramn2(ii) ) num2str( paramn3(ii) )  ...
    '_ls' num2str( paraml1(ii) ) num2str( paraml2(ii) )...
    '_bc' num2str( parambc(ii), '%.2f' ) '_vD' num2str( paramvD(ii), '%.3g' ) ...
    noiseStr{ii}, ...
    '_diffNl' num2str( nlDiff, '%.2f' ), ...
    '_IC' initStr '_SM' num2str( paramSM(ii), '%d'  ) ...
    '_t' num2str( trial,'%.2d' ) '.' num2str( paramRun(ii), '%.2d' ) '.mat' ];
  fprintf('\nStarting %s \n', filename);
  [denRecObjTemp] = ddftMain( filename, paramvec, systemObj, particleObj,...
    runObj, timeObj, rhoInitObj, flags, ...
    interactParam, extParam, noiseParam);
  denRecObj{ii} = denRecObjTemp;
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
