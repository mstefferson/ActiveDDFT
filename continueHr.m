function [denRecObj] = continueHr()
% some global
global runSave
try
  dateTime =  datestr(now);
  fprintf('Starting RunHardRod: %s\n', dateTime);
  ticID = tic;
  % Add Subroutine path
  currentDir = pwd;
  addpath( genpath( [currentDir '/src'] ) );
  evolvedSucess = 0;
  % Record how long things take
  tMainID  = tic;
  % find parameters/runfiles. Make sure there is only one
  runlist = dir('./run_*');
  paramslist = dir('./params_*');
  if length( runlist ) > 1 || length( paramslist ) > 1
    error('too many runs to finish!')
  elseif isempty( runlist ) || isempty( paramslist )
    error('no runs to finish!')
  end
  % get names
  filename = runlist.name( 5:end );
  saveNameRun   = ['run_' filename];
  saveNameRhoFinal   = ['rhoFinal_' filename];
  saveNameParams = ['params_' filename];
  % Create a file that holds warning print statements
  locString = sprintf('Loc_%s.txt', filename(1:end-4));
  lfid      = fopen(locString,'a+');    % a+ allows to append data
  fprintf(lfid,'Starting main, current code\n');
  % load params, matfile runs
  load( saveNameParams )
  paramSave = matfile(saveNameParams,'Writable',true);
  runSave = matfile(runlist.name,'Writable',true);
  rhoFinalSave = matfile(saveNameRhoFinal,'Writable',true);
  % get directory name
  dirName  = denRecObj.dirName;
  dirId  = denRecObj.dirId;
  if flags.MakeOP
    saveNameOP   = ['op_' filename];
    opSave = matfile(saveNameOP,'Writable',true);
  end
  if exist(dirName,'dir') == 0
    mkdir(dirName);
  end
  % Make remaining objects
  [timeObj, gridObj, diffObj, interObj, noise, polarDrive, ...
    dRhoFlux, densityDepDiff, runTime] = buildClassesAndStructures( ...
    systemObj, particleObj, flags, timeObj, lfid );
  % get run time
  timeObjCont = timeObj;
  % Fix
  ntCompleted = size( runSave.Den_rec, 4 );
  % be careful if ntComplete = 2 b/c it could be zeros!!!
  if ntCompleted == 2; ntCompleted = 1; end
  timeObjCont.recStartInd = ntCompleted + 1;
  nRecNew = timeObj.N_rec - ntCompleted;
  timeObjCont.N_rec = nRecNew;
  timeObjCont.N_time = timeObjCont.N_rec * timeObj.N_dtRec;
  % grab rho
  rho = runSave.Den_rec(:,:,:,ntCompleted);
  % print some things
  myStr = [ 'Ran ' num2str( ntCompleted, '%d') ' chunks. ' ...
    'Want ' num2str( timeObj.N_rec, '%d') ' total. Running ' ...
    num2str( timeObjCont.N_time, '%d') ' more time steps starting recording at ' ...
    num2str( timeObjCont.recStartInd, '%d')];
  fprintf('%s\n', myStr )
  % Create a file that holds warning print statements
  locString = sprintf('Loc_%s.txt', filename(1:end-4));
  lfid = fopen(locString,'a+');    % a+ allows to append data
  % Run the main code
  tBodyID      = tic;
  if ntCompleted < timeObj.N_rec
    fprintf('Continue run over time\n');
    if flags.DiagLop == 1
      [denRecObj, rho]  = denEvolverFTDiagOp(...
        rho, systemObj, timeObjCont, gridObj, diffObj, interObj, ...
        polarDrive, noise, dRhoFlux, densityDepDiff, flags, lfid);
    else
      [denRecObj, rho]  = denEvolverFT(...
        rho, systemObj, timeObjCont, gridObj, diffObj, interObj, ...
        polarDrive, noise, dRhoFlux, densityDepDiff, flags, lfid);
    end
    % Save it
    if flags.SaveMe
      paramSave.denRecObj = denRecObj;
      runSave.denRecObj   = denRecObj;
      rhoFinalSave.rho   = rho;
    end
  else
    fprintf('Run over time finished. Just need ops\n');
    denRecObj = paramSave.denRecObj;
  end
  bodyRunTime  = toc(tBodyID);
  evolvedSucess = 1;
  denRecObj.dirName = dirName;
  if flags.Verbose
    fprintf('Ran Main Body t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, bodyRunTime);
  end
  fprintf(lfid,'Body Run Time = %f\n\n', bodyRunTime);
  runTime.body = bodyRunTime;
  % Make order params if you want
  if flags.MakeOP  == 1
    tOpID           = tic ;
    opSave = buildOPs(saveNameOP, systemObj, gridObj, denRecObj, ...
      timeObj, rhoInit, runSave);
    opRunTime = toc(tOpID);
    if flags.Verbose
      fprintf('Made OP object t%d_%d: %.3g \n', ...
        runObj.trialID, runObj.runID, opRunTime);
    end
    fprintf(lfid,'OrderParam Run time = %f\n', opRunTime);
    runTime.op = opRunTime;
  end %if OP
  % Move saved things
  if flags.SaveMe
    paramSave.runTime = runTime;
    movefile(saveNameRun,dirName);
    movefile( saveNameRhoFinal,dirName);
    movefile( saveNameParams,dirName);
    if flags.MakeOP == 1
      movefile( saveNameOP,dirName);
    end
  end
  % run analysis
  if flags.movie.analysis == 1
    movieHardRod(dirId)
  end
  % Save how long everything took
  totRunTime = toc(tMainID);
  if flags.Verbose
    fprintf('Run Finished t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, totRunTime);
  end
  fprintf(lfid,'Total Run time = %f\n', totRunTime);
  runTime.tot = totRunTime;
  % move file to indicate that you are done
  movefile('ParamsRunning.mat', 'ParamsFinished.mat');
catch err %Catch errors
  % write the error to file and to screen
  fprintf('%s', err.getReport('extended')) ;
  runSave.err = err;
end %End try and catch
% close files
if exist('lfid','var')
  fclose(lfid);
end
if exist('flags','var')
  if flags.Verbose
    fprintf('Leaving Main for t%d.%d\n', ...
      runObj.trialID, runObj.runID);
  end
end
runTime = toc(ticID);
dateTime =  datestr(now);
runHr = floor( runTime / 3600); runTime = runTime - runHr*3600;
runMin = floor( runTime / 60);  runTime = runTime - runMin*60;
runSec = floor(runTime);
fprintf('RunTime: %.2d:%.2d:%.2d (hr:min:sec)\n', runHr, runMin,runSec);
fprintf('Finished RunHardRod: %s\n', dateTime);
end
