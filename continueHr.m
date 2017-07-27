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
  % load params, matfile runs
  load( saveNameParams )
  paramSave = matfile(saveNameParams,'Writable',true);
  runSave = matfile(runlist.name,'Writable',true);
  rhoFinalSave = matfile(saveNameRhoFinal,'Writable',true);
  % get directory name
  dirName  = filename(1:end-4) ;
  if flags.MakeOP == 0
    dirName  = ['./runfiles/' dirName ];
  else
    saveNameOP   = ['op_' filename];
    if flags.MakeMovies == 1
      dirName  = ['./analyzedfiles/' dirName ];
    else
      dirName  = ['./runOPfiles/' dirName ];
    end
    opSave = matfile(saveNameOP,'Writable',true);
  end
  if exist(dirName,'dir') == 0
    mkdir(dirName);
  end
  %dirName = denRecObj.dirName;
  % get run time
  timeObjCont = timeObj;
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
  fprintf('Ran %d chunks. Want %d total. Running %d more time steps starting recording at %d\n', ...
    ntCompleted, timeObj.N_rec, timeObjCont.N_time, timeObjCont.recStartInd )
  % Create a file that holds warning print statements
  locString = sprintf('Loc_%s.txt', filename(1:end-4));
  lfid      = fopen(locString,'a+');    % a+ allows to append data
  % rerun some unsaved things
  % grid
  tGridID = tic;
  [gridObj] = GridMakerPBCxk(systemObj.n1,systemObj.n2,systemObj.n3,...
    systemObj.l1,systemObj.l2,systemObj.l3);
  gridrunTime = toc(tGridID);
  runTime.grid = gridrunTime;
  % diff
  tDiffID = tic;
  [diffObj] =  DiffMobCoupCoeffCalc( systemObj.tmp,...
    particleObj.mob,particleObj.mobPar,particleObj.mobPerp,particleObj.mobRot,...
    gridObj.k1, gridObj.k2, gridObj.k3, ...
    gridObj.k1rep2, gridObj.k2rep2,particleObj.vD);
  diffRunTime = toc(tDiffID);
  runTime.diff = diffRunTime;
  [interObj] =  interObjMaker( particleObj, systemObj, gridObj );
  % Run the main code
  tBodyID      = tic;
  if flags.DiagLop == 1
    [denRecObj, rho]  = denEvolverFTDiagOp(...
      rho, systemObj, particleObj, timeObjCont, gridObj, ...
      diffObj, interObj, flags, lfid);
  else
    [denRecObj, rho]  = denEvolverFT(...
      rho, systemObj, particleObj, timeObjCont, gridObj, ...
      diffObj, interObj, flags, lfid);
  end
  bodyRunTime  = toc(tBodyID);
  evolvedSucess = 1;
  denRecObj.dirName = dirName;
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
    tOpID           = tic ;
    % Save
    if systemObj.n3 > 1
      % Commonly used trig functions
      [~,~,phi3D] = meshgrid(gridObj.x2,gridObj.x1,gridObj.x3);
      cosPhi3d = cos(phi3D);
      sinPhi3d = sin(phi3D);
      cos2Phi3d = cosPhi3d .^ 2;
      sin2Phi3d = sinPhi3d .^ 2;
      cossinPhi3d = cosPhi3d .* sinPhi3d;
    else
      cosPhi3d=0;sinPhi3d=0;cos2Phi3d=0;sin2Phi3d=0;cossinPhi3d=0;
    end
    % Build time rec vector
    if denRecObj.DidIBreak == 0
      totRec = length( denRecObj.TimeRecVec );
      opTimeRecVec = denRecObj.TimeRecVec ;
      opSave.OpTimeRecVec = opTimeRecVec;
    else %Don't incldue the blowed up denesity for movies. They don't like it.
      totRec = length( denRecObj.TimeRecVec ) - 1;
      opTimeRecVec = denRecObj.TimeRecVec(1:end-1) ;
      opSave.OpTimeRecVec = opTimeRecVec;
    end
    % Set up saving
    opSave.C_rec    = zeros(systemObj.n1, systemObj.n2, 2);
    if systemObj.n3 > 1
      % Distribution slice
      holdX = systemObj.n1 /2 + 1; % spatial pos placeholders
      holdY = systemObj.n2 /2 + 1; % spatial pos placeholders
      opSave.distSlice_rec = reshape( ...
        runSave.Den_rec(holdX, holdY, : , 1:length(opTimeRecVec)),...
        [systemObj.n3 length(opTimeRecVec)] );
      opSave.POP_rec  = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.POPx_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.POPy_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.NOP_rec  = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.NOPx_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.NOPy_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.aveC_rec = zeros(1,2);
      opSave.aveP_rec = zeros(1,2);
      opSave.aveN_rec = zeros(1,2);
    end
    % Break it into chunks
    numChunks = timeObj.N_chunks;
    sizeChunk = floor( totRec/ numChunks );
    if sizeChunk > 0
      numChunks = ceil( totRec/ sizeChunk);
    else
      numChunks = 1;
    end
    for i = 1:numChunks
      if i ~= numChunks
        currInd =  (i-1) * sizeChunk + 1: i * sizeChunk;
      else
        if numChunks == 1
          currInd = 1:totRec;
        else
          currInd = (i-1) * sizeChunk:totRec;
        end
      end
      % Make the records
      [OPObjTemp] = CPNrecMaker(systemObj.n1,systemObj.n2,...
        opTimeRecVec(currInd), runSave.Den_rec(:,:,:,currInd) ,...
        gridObj.x3,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d );
      % Save it
      opSave.C_rec(:,:,currInd) = OPObjTemp.C_rec;
      if systemObj.n3 > 1
        opSave.POP_rec(:,:,currInd) = OPObjTemp.POP_rec;
        opSave.POPx_rec(:,:,currInd) = OPObjTemp.POPx_rec;
        opSave.POPy_rec(:,:,currInd) = OPObjTemp.POPy_rec;
        opSave.NOP_rec(:,:,currInd) = OPObjTemp.NOP_rec;
        opSave.NOPx_rec(:,:,currInd) = OPObjTemp.NOPx_rec;
        opSave.NOPy_rec(:,:,currInd) = OPObjTemp.NOPy_rec;
        opSave.aveC_rec(1,currInd) = OPObjTemp.aveC_rec;
        opSave.aveP_rec(1,currInd) = OPObjTemp.aveP_rec;
        opSave.aveN_rec(1,currInd) = OPObjTemp.aveN_rec;
   end
    end % loop over chunks
    % fix size
    opSave.C_rec    = opSave.C_rec(:,:,1:totRec);
    opSave.aveC_rec = opSave.aveC_rec(1,1:totRec);
    if systemObj.n3 > 1
      opSave.POP_rec  = opSave.POP_rec(:,:,1:totRec);
      opSave.POPx_rec = opSave.POPx_rec(:,:,1:totRec);
      opSave.POPy_rec = opSave.POPy_rec(:,:,1:totRec);
      opSave.NOP_rec  = opSave.NOP_rec(:,:,1:totRec);
      opSave.NOPx_rec = opSave.NOPx_rec(:,:,1:totRec);
      opSave.NOPy_rec = opSave.NOPy_rec(:,:,1:totRec);
      opSave.aveP_rec = opSave.aveP_rec(1,1:totRec);
      opSave.aveN_rec = opSave.aveN_rec(1,1:totRec);
    end
    if flags.MakeMovies
      OPobj.OpTimeRecVec = opTimeRecVec;
      OPobj.C_rec    = opSave.C_rec;
      if systemObj.n3 > 1
        OPobj.POP_rec  = opSave.POP_rec;
        OPobj.POPx_rec = opSave.POPx_rec;
        OPobj.POPy_rec = opSave.POPy_rec;
        OPobj.NOP_rec  = opSave.NOP_rec;
        OPobj.NOPx_rec = opSave.NOPx_rec;
        OPobj.NOPy_rec = opSave.NOPy_rec;
      end
    end
    % Now do it for steady state sol
    if systemObj.n3 > 1
      [~,~,phi3D] = meshgrid(1,1,gridObj.x3);
      cosPhi3d = cos(phi3D);
      sinPhi3d = sin(phi3D);
      cos2Phi3d = cosPhi3d .^ 2;
      sin2Phi3d = sinPhi3d .^ 2;
      cossinPhi3d = cosPhi3d .* sinPhi3d;
      % Calc CPN
      [~,~,~,~,opSave.NOPeq,~,~] = ...
        OpCPNCalc(1, 1, reshape( rhoInit.feq, [1,1,systemObj.n3] ), ...
        gridObj.x3,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);
      % Save if movies
      if flags.MakeMovies
        OPobj.NOPeq = opSave.NOPeq;
      end
    end
    opRunTime = toc(tOpID);
    if flags.Verbose
      fprintf('Made OP object t%d_%d: %.3g \n', ...
        runObj.trialID, runObj.runID, opRunTime);
    end
    fprintf(lfid,'OrderParam Run time = %f\n', opRunTime);
    runTime.op = opRunTime;
    % Make movies
    if flags.MakeMovies == 1
      movieSuccess = 0;
      % Make matlab movies
      tMovID       = tic;
      % for n3 = 1, just concentration
      if systemObj.n3 == 1
        % Save Name
        movStr = sprintf('Cmov_bc%.2f_vD%.1f_%.2d_%.2d.avi',...
          systemObj.bc,particleObj.vD,runObj.trialID, runObj.runID);
        % Run function
        CMovieMakerAvi(movStr,...
          gridObj.x1,gridObj.x2,particleObj.b .* OPobj.C_rec,...
          OPobj.OpTimeRecVec);
      else
        % Save Name
        movStr = sprintf('OPmov_bc%.2f_vD%.1f_%.2d_%.2d.avi',...
          systemObj.bc,particleObj.vD,runObj.trialID, runObj.runID);
        % Run function
        OPMovieMakerTgtherDirAvi(movStr,...
          gridObj.x1,gridObj.x2,gridObj.x3,OPobj,...
          OPobj.distSlice_rec,OPobj.OpTimeRecVec);
      end
      movieSuccess = 1;
      % Move it
      movRunTime   = toc(tMovID);
      if flags.Verbose
        fprintf('Made movies t%d_%d: %.3g \n', ...
          runObj.trialID, runObj.runID, movRunTime);
      end
      fprintf(lfid,'Make Mov Run Time = %f\n',  movRunTime);
      runTime.mov = movRunTime;
      % Make amplitude plot
      k1ind0 = gridObj.k1ind0;
      k2ind0 = gridObj.k2ind0;
      k3ind0 = gridObj.k3ind0;
      nRec = length( denRecObj.TimeRecVec);
      % Ft amps
      if systemObj.n3 > 1
        totModes   = 12;
      else
        totModes = 4;
      end
      FTind2plot = zeros( totModes , 3 );
      FTmat2plot = zeros( totModes , nRec );
      FTind2plot(1,:) = [k1ind0     k2ind0     k3ind0 ];
      FTind2plot(2,:) = [k1ind0 + 1 k2ind0     k3ind0 ];
      FTind2plot(3,:) = [k1ind0     k2ind0 + 1 k3ind0 ];
      FTind2plot(4,:) = [k1ind0 + 1 k2ind0 + 1 k3ind0 ];
      if systemObj.n3 > 1
        FTind2plot(5,:) = [k1ind0     k2ind0     k3ind0 + 1];
        FTind2plot(6,:) = [k1ind0 + 1 k2ind0     k3ind0 + 1];
        FTind2plot(7,:) = [k1ind0     k2ind0 + 1 k3ind0 + 1];
        FTind2plot(8,:) = [k1ind0 + 1 k2ind0 + 1 k3ind0 + 1];
        FTind2plot(9,:) = [k1ind0     k2ind0     k3ind0 + 2];
        FTind2plot(10,:) = [k1ind0 + 1 k2ind0     k3ind0 + 2];
        FTind2plot(11,:) = [k1ind0     k2ind0 + 1 k3ind0 + 2];
        FTind2plot(12,:) = [k1ind0 + 1 k2ind0 + 1 k3ind0 + 2];
      end
      % Scale by N so it's N independent
      for i = 1:totModes
        FTmat2plot(i,:) =  1 / (systemObj.n1 * systemObj.n2 * systemObj.n3) .* ...
          reshape(runSave.DenFT_rec( FTind2plot(i,1), FTind2plot(i,2), FTind2plot(i,3),1:nRec ),...
          [ 1, nRec ]  );
      end
      % Plot Amplitudes
      ampPlotterFT(FTmat2plot, FTind2plot, ...
        denRecObj.TimeRecVec(1:nRec), k1ind0, k2ind0, k3ind0, timeObj.t_tot);
      % Save it
      figtl = sprintf('AmpFT.fig');
      % savefig doesn't like decimals so save it and rename it.
      savefig(gcf,figtl)
      figtl2 = sprintf('AmpFT_bc%.2f_vD%.0f_%.2d_%.2d',...
        systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
      movefile(figtl,[figtl2 '.fig'])
      saveas(gcf, [figtl2 '.jpg'],'jpg')
      % Plot final slices of final order parameters
      if systemObj.n3 > 1
        sliceSaveTag = sprintf('SOP_bc%.2f_vD%.0f_%.2d_%.2d',...
          systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
        if denRecObj.DidIBreak == 0
          sliceOPplot( OPobj.C_rec(:,:,end), OPobj.POP_rec(:,:,end),...
            OPobj.NOP_rec(:,:,end), systemObj, ...
            gridObj, rhoFinalSave.rhoFinal, sliceSaveTag )
        else
          stepsNb = length(OPobj.OpTimeRecVec);
          sliceOPplot( OPobj.C_rec(:,:,end), OPobj.POP_rec(:,:,end),...
            OPobj.NOP_rec(:,:,end), systemObj, ...
            gridObj, runSave.Den_rec( :,:,:, stepsNb), sliceSaveTag );
        end
        % Move it
        movefile([sliceSaveTag '*'], dirName);
        % Plot max order parameters vs time
        maxSaveTag = sprintf('MaxOP_bc%.2f_vD%.0f_%.2d_%.2d',...
          systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
        plotMaxOPvsTime( OPobj.C_rec, OPobj.POP_rec, OPobj.NOP_rec, ...
          particleObj.b, OPobj.OpTimeRecVec, maxSaveTag );
      else
        % Plot max order parameters vs time
        maxSaveTag = sprintf('MaxC_bc%.2f_vD%.0f_%.2d_%.2d',...
          systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
        plotMaxCvsTime( OPobj.C_rec, particleObj.b, OPobj.OpTimeRecVec, maxSaveTag );
      end
      % Move everything else
      movefile([figtl2 '*'], dirName);
      movefile([movStr '*'], dirName);
      movefile([maxSaveTag '*'], dirName);
    end % End if movies
  end % if OP
  % Save how long everything took
  totRunTime = toc(tMainID);
  if flags.Verbose
    fprintf('Run Finished t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, totRunTime);
  end
  fprintf(lfid,'Total Run time = %f\n', totRunTime);
  runTime.tot = totRunTime;
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
catch err %Catch errors
  % write the error to file and to screen
  fprintf('%s', err.getReport('extended')) ;
  runSave.err = err;
  % Movies can have issues to box size. If they do, just move files
  % to ./runOPfiles
  % Move saved things
  if evolvedSucess == 1 && flags.SaveMe == 1
    if flags.MakeMovies == 1 && movieSuccess == 0
      fprintf('Movies failed\n');
      if exist(movStr,'file'); delete(movStr); end
      rmdir(dirName);
      dirName = [ './runOPfiles' filename(1:end-4) ];
      denRecObj.dirName = dirName;
    end
    movefile(saveNameRun,dirName);
    movefile( saveNameRhoFinal,dirName);
    movefile( saveNameParams,dirName);
    if flags.MakeOP == 1
      movefile( saveNameOP,dirName);
    end
  end
end %End try and catch
% close files
fclose(lfid);
delete(locString);
if flags.Verbose
  fprintf('Leaving Main for t%d.%d\n', ...
    runObj.trialID, runObj.runID);
end
runTime = toc(ticID);
dateTime =  datestr(now);
movefile('ParamsRunning.mat', 'ParamsFinished.mat');
runHr = floor( runTime / 3600); runTime = runTime - runHr*3600;
runMin = floor( runTime / 60);  runTime = runTime - runMin*60;
runSec = floor(runTime);
fprintf('RunTime: %.2d:%.2d:%.2d (hr:min:sec)\n', runHr, runMin,runSec);
fprintf('Finished RunHardRod: %s\n', dateTime);
end
