
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT

function  [ denRecObj ] = ...
  ddftMain( filename, paramVec, systemObj, particleObj, runObj, timeObj, ...
  rhoInit, flags, interactLrV, externalV, noise )
% use latex for plots
set(0,'defaulttextinterpreter','latex')
% some global
global runSave
try
  % Create a file that holds warning print statements
  locString = sprintf('Loc_%s.txt', filename(1:end-4));
  lfid      = fopen(locString,'a+');    % a+ allows to append data
  fprintf(lfid,'Starting main, current code\n');
  % Set up denRecObj just in case code doesn't finish
  denRecObj.didIrun = 0;
  % Move parameter vector to obj
  systemObj.n1 = paramVec(1);
  systemObj.n2 = paramVec(2);
  systemObj.n3 = paramVec(3);
  systemObj.l1 = paramVec(4);
  systemObj.l2 = paramVec(5);
  particleObj.vD = paramVec(6);
  systemObj.bc = paramVec(7);
  flags.StepMeth = paramVec(8);
  runObj.runID = paramVec(9);
  particleObj.interactLrV = interactLrV{paramVec(10)};
  particleObj.externalV = externalV{paramVec(11)};
  systemObj.noise = noise{paramVec(12)};
  % set some other parameters
  systemObj.c = systemObj.bc ./ particleObj.b;
  systemObj.numPart  = systemObj.c * systemObj.l1 * systemObj.l2; % number of particles
  systemObj.lBox = [systemObj.l1 systemObj.l2];
  % Set-up save paths, file names, and matfile
  if flags.SaveMe
    saveNameRun   = ['run_' filename];
    saveNameRhoFinal   = ['rhoFinal_' filename];
    saveNameParams = ['params_' filename];
    dirName  = filename(1:end-4) ;
    if flags.MakeOP == 0
      dirPath = 'runfiles/';
    else
      saveNameOP   = ['op_' filename];
      dirPath = 'runOPfiles/';
      opSave = matfile(saveNameOP,'Writable',true);
    end
    % dont overwrite dirs
    if exist( [dirPath dirName],'dir') == 0
      dirId = dirName;
    else
      dirId = [dirName '_' datestr(now,'yyyymmdd_HH.MM.SS') ];
    end
    dirName = [dirPath dirId];
    mkdir( dirName );
    paramSave = matfile(saveNameParams,'Writable',true);
    rhoFinalSave = matfile(saveNameRhoFinal,'Writable',true);
    runSave = matfile(saveNameRun,'Writable',true);
    denRecObj.dirName = dirName; % Just in case
    denRecObj.dirId = dirId; % Just in case
  else
    dirName = pwd;
  end
  % verbose
  if flags.Verbose
    if flags.AnisoDiff
      fprintf('Running Anisotropic Diffusion\n');
    else
      fprintf('Running Isotropic Diffusion\n');
    end
  end
  % Record how long things take
  tMainID  = tic;
  % Make remaining objects
  [timeObj, gridObj, diffObj, interObj, noise, polarDrive, dRhoFlux, runTime] ...
    = buildClassesAndStructures( systemObj, particleObj, ...
    flags, timeObj, lfid );
  % Build initial density
  %Initialze density
  tIntDenID = tic;
  [rho] = MakeConc(systemObj,particleObj,rhoInit,gridObj);
  intDenRunTime = toc(tIntDenID);
  if flags.Verbose
    fprintf('Made initial density t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, intDenRunTime);
  end
  fprintf(lfid,'Made initial density: %.3g\n', intDenRunTime);
  runTime.intDen = intDenRunTime;  % Save everything before running body of code
  flags.dRhoCalc = interObj.anyInter || polarDrive.Flag || ...
    noise.Flag;
  if flags.SaveMe
    runSave.flags    = flags;
    runSave.runObj    = runObj;
    runSave.systemObj = systemObj;
    runSave.particleObj = particleObj;
    runSave.timeObj  = timeObj;
    runSave.rhoInit  = rhoInit;
    % Clean up gridobj before saving
    fields2del = {'k1rep2','k2rep2'};
    gridTemp = rmfield(gridObj,fields2del);
    runSave.gridObj  = gridTemp;
    runSave.Den_rec = zeros(systemObj.n1,systemObj.n2,systemObj.n3,2);
    runSave.DenFT_rec = complex( ...
      zeros(systemObj.n1,systemObj.n2,systemObj.n3,2), 0 );
    runSave.Den_rec(:,:,:,1) = rho;
    runSave.DenFT_rec(:,:,:,1) = fftshift(fftn(rho));
    runSave.denRecObj   = denRecObj;
    runSave.numSavedRhos = 1;
    % Save params now
    paramSave.flags = flags;
    paramSave.particleObj = particleObj;
    paramSave.rhoInit = rhoInit;
    paramSave.systemObj = systemObj;
    paramSave.runObj = runObj;
    paramSave.timeObj = timeObj;
    paramSave.denRecObj = denRecObj;
  end
  % Run the main code
  tBodyID      = tic;
  if flags.DiagLop == 1
    [denRecObj, rho]  = denEvolverFTDiagOp( ...
      rho, systemObj, timeObj, gridObj, diffObj, interObj, ...
      polarDrive, noise, dRhoFlux, flags, lfid);
  else
    [denRecObj, rho]  = denEvolverFT( ...
      rho, systemObj, timeObj, gridObj, diffObj, interObj, ...
      polarDrive, noise, dRhoFlux, flags, lfid);
  end
  bodyRunTime  = toc(tBodyID);
  % Save it
  if flags.SaveMe
    paramSave.denRecObj = denRecObj;
    runSave.denRecObj   = denRecObj;
    rhoFinalSave.rho   = rho;
  end
  if flags.Verbose
    fprintf('Ran Main Body t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, bodyRunTime);
  end
  fprintf(lfid,'Body Run Time = %f\n\n', bodyRunTime);
  runTime.body = bodyRunTime;
  % Run movies if you want
  if flags.MakeOP  == 1
    tOpID = tic ;
    opSave = buildOPs(saveNameOP, systemObj, gridObj, denRecObj, ...
      timeObj, rhoInit, runSave);
    opRunTime = toc(tOpID);
    if flags.Verbose
      fprintf('Made OP object t%d_%d: %.3g \n', ...
        runObj.trialID, runObj.runID, opRunTime);
    end
    fprintf(lfid,'OrderParam Run time = %f\n', opRunTime);
    runTime.op = opRunTime;
  end % if OP
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
  end  % Save how long everything took
  totRunTime = toc(tMainID);
  if flags.Verbose
    fprintf('Run Finished t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, totRunTime);
  end
  fprintf(lfid,'Total Run time = %f\n', totRunTime);
  runTime.tot = totRunTime;
catch err %Catch errors
  % write the error to file and to screen
  fprintf('%s\n', err.getReport('extended')) ;
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
end % End ddftMain
